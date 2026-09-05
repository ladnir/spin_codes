#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <intrin.h>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
#ifndef INDEPENDENT_VERIFIER
#define INDEPENDENT_VERIFIER 0
#endif
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned Nodes = 32772;
    constexpr unsigned Coordinates = 64;
    constexpr std::array<unsigned, 7> Degrees{1, 2, 4, 9, 10, 18, 20};
    constexpr std::array<std::uint64_t, 7> Factors{
        0x3ULL, 0x7ULL, 0x13ULL, 0x373ULL, 0x519ULL, 0x7C9C3ULL, 0x1E1FFFULL};
    constexpr std::array<std::uint32_t, 7> Orders{
        1, 3, 15, 511, 1023, 37449, 1048575};
#if INDEPENDENT_VERIFIER
    constexpr std::array<unsigned, 7> CanonicalOrder{0, 1, 2, 3, 4, 5, 6};
#else
    constexpr std::array<unsigned, 7> CanonicalOrder{6, 5, 4, 3, 2, 1, 0};
#endif
    constexpr std::uint64_t ExpectedRecordCount =
        Coordinates * (Coordinates - 1) / 2 +
        std::uint64_t{Nodes - 1} * Coordinates * Coordinates;

    struct U128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;

        U128() = default;
        U128(std::uint64_t value) : low(value) {}

        U128& operator+=(const U128& other)
        {
            const std::uint64_t old = low;
            low += other.low;
            high += other.high + (low < old);
            return *this;
        }
    };

    U128 operator+(U128 left, const U128& right)
    {
        left += right;
        return left;
    }

    bool operator==(const U128& left, const U128& right)
    {
        return left.low == right.low && left.high == right.high;
    }

    U128 multiply(std::uint64_t left, std::uint64_t right)
    {
        U128 result;
        result.low = _umul128(left, right, &result.high);
        return result;
    }

    U128 operator*(const U128& left, std::uint64_t right)
    {
        std::uint64_t carry = 0;
        U128 result;
        result.low = _umul128(left.low, right, &carry);
        const std::uint64_t upper = left.high * right;
        result.high = carry + upper;
        if (result.high < carry)
            throw std::overflow_error("dual A4: 128-bit multiplication overflow");
        return result;
    }

    U128 operator*(std::uint64_t left, const U128& right)
    {
        return right * left;
    }

    unsigned remainder(const U128& value, unsigned divisor)
    {
        const unsigned highRemainder = static_cast<unsigned>(value.high % divisor);
        std::uint64_t ignored = 0;
        (void)_udiv128(highRemainder, value.low, divisor, &ignored);
        return static_cast<unsigned>(ignored);
    }

    U128 divide(U128 value, unsigned divisor, unsigned& finalRemainder)
    {
        U128 quotient;
        quotient.high = value.high / divisor;
        const std::uint64_t highRemainder = value.high % divisor;
        std::uint64_t lowRemainder = 0;
        quotient.low = _udiv128(highRemainder, value.low, divisor, &lowRemainder);
        finalRemainder = static_cast<unsigned>(lowRemainder);
        return quotient;
    }

    std::string decimal(U128 value)
    {
        if (value.high == 0 && value.low == 0)
            return "0";
        std::string result;
        while (value.high || value.low)
        {
            unsigned digit = 0;
            value = divide(value, 10, digit);
            result.push_back(static_cast<char>('0' + digit));
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::uint64_t carrylessMultiplyLow(std::uint64_t left, std::uint64_t right)
    {
        std::uint64_t result = 0;
        while (right)
        {
            const unsigned bit = std::countr_zero(right);
            result ^= left << bit;
            right &= right - 1;
        }
        return result;
    }

    std::uint64_t generatorInverse()
    {
        std::uint64_t inverse = 1;
        for (unsigned bit = 1; bit < 64; ++bit)
        {
            if ((carrylessMultiplyLow(Generator, inverse) >> bit) & 1)
                inverse |= std::uint64_t{1} << bit;
        }
        if (carrylessMultiplyLow(Generator, inverse) != 1)
            throw std::runtime_error("dual A4: generator inverse failed");
        return inverse;
    }

    struct Word128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };

    Word128 encodeBch(std::uint64_t message)
    {
        Word128 result{};
        while (message)
        {
            const unsigned bit = std::countr_zero(message);
            result.low ^= Generator << bit;
            result.high ^=
                (bit == 0 ? 0 : Generator >> (64 - bit)) |
                (std::uint64_t{1} << 63);
            message &= message - 1;
        }
        return result;
    }

    std::array<std::uint64_t, 64> systematicRightColumns()
    {
        const std::uint64_t inverse = generatorInverse();
        std::array<std::uint64_t, 64> result{};
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            const std::uint64_t message = carrylessMultiplyLow(
                std::uint64_t{1} << bit, inverse);
            const Word128 word = encodeBch(message);
            if (word.low != (std::uint64_t{1} << bit))
                throw std::runtime_error("dual A4: systematic BCH reconstruction failed");
            result[bit] = word.high;
        }
        return result;
    }

    std::uint64_t accumulate(std::uint64_t value)
    {
        value ^= value << 1;
        value ^= value << 2;
        value ^= value << 4;
        value ^= value << 8;
        value ^= value << 16;
        value ^= value << 32;
        return value;
    }

    class ByteMap
    {
    public:
        explicit ByteMap(const std::array<std::uint64_t, 64>& columns)
        {
            for (unsigned byteIndex = 0; byteIndex < 8; ++byteIndex)
            {
                for (unsigned byte = 0; byte < 256; ++byte)
                {
                    for (unsigned bit = 0; bit < 8; ++bit)
                    {
                        if ((byte >> bit) & 1)
                            mTables[byteIndex][byte] ^= columns[8 * byteIndex + bit];
                    }
                }
            }
        }

        std::uint64_t operator()(std::uint64_t value) const
        {
            return
                mTables[0][value & 0xff] ^
                mTables[1][(value >> 8) & 0xff] ^
                mTables[2][(value >> 16) & 0xff] ^
                mTables[3][(value >> 24) & 0xff] ^
                mTables[4][(value >> 32) & 0xff] ^
                mTables[5][(value >> 40) & 0xff] ^
                mTables[6][(value >> 48) & 0xff] ^
                mTables[7][value >> 56];
        }

    private:
        std::array<std::array<std::uint64_t, 256>, 8> mTables{};
    };

    std::array<std::uint64_t, 64> transposeColumns(
        const std::array<std::uint64_t, 64>& columns)
    {
        std::array<std::uint64_t, 64> result{};
        for (unsigned output = 0; output < 64; ++output)
        {
            for (unsigned input = 0; input < 64; ++input)
                result[output] |= ((columns[input] >> output) & 1) << input;
        }
        return result;
    }

    std::uint64_t applyPolynomial(
        const ByteMap& step,
        std::uint64_t polynomial,
        std::uint64_t value)
    {
        std::uint64_t result = 0;
        std::uint64_t current = value;
        while (polynomial)
        {
            if (polynomial & 1)
                result ^= current;
            polynomial >>= 1;
            current = step(current);
        }
        return result;
    }

    std::vector<std::uint64_t> kernelBasis(
        const ByteMap& step,
        std::uint64_t polynomial,
        unsigned expectedDimension)
    {
        std::array<std::uint64_t, 64> pivots{};
        std::array<std::uint64_t, 64> representations{};
        std::vector<std::uint64_t> result;
        for (unsigned column = 0; column < 64; ++column)
        {
            std::uint64_t value = applyPolynomial(
                step, polynomial, std::uint64_t{1} << column);
            std::uint64_t representation = std::uint64_t{1} << column;
            while (value)
            {
                const unsigned pivot = std::bit_width(value) - 1;
                if (pivots[pivot])
                {
                    value ^= pivots[pivot];
                    representation ^= representations[pivot];
                }
                else
                {
                    pivots[pivot] = value;
                    representations[pivot] = representation;
                    break;
                }
            }
            if (value == 0)
                result.push_back(representation);
        }
        if (result.size() != expectedDimension)
            throw std::runtime_error("dual A4: component kernel dimension mismatch");
        return result;
    }

    class CoordinateMap
    {
    public:
        explicit CoordinateMap(const std::vector<std::uint64_t>& basis)
        {
            if (basis.size() != 64)
                throw std::runtime_error("dual A4: coordinate basis is not full");
            std::array<std::uint64_t, 64> pivots{};
            std::array<std::uint64_t, 64> representations{};
            for (unsigned coordinate = 0; coordinate < 64; ++coordinate)
            {
                std::uint64_t value = basis[coordinate];
                std::uint64_t representation = std::uint64_t{1} << coordinate;
                while (value)
                {
                    const unsigned pivot = std::bit_width(value) - 1;
                    if (pivots[pivot])
                    {
                        value ^= pivots[pivot];
                        representation ^= representations[pivot];
                    }
                    else
                    {
                        pivots[pivot] = value;
                        representations[pivot] = representation;
                        break;
                    }
                }
                if (value == 0)
                    throw std::runtime_error("dual A4: dependent coordinate basis");
            }
            std::array<std::uint64_t, 64> columns{};
            for (unsigned bit = 0; bit < 64; ++bit)
            {
                std::uint64_t value = std::uint64_t{1} << bit;
                std::uint64_t result = 0;
                while (value)
                {
                    const unsigned pivot = std::bit_width(value) - 1;
                    if (!pivots[pivot])
                        throw std::runtime_error("dual A4: coordinate solve failed");
                    value ^= pivots[pivot];
                    result ^= representations[pivot];
                }
                columns[bit] = result;
            }
            mMap = new ByteMap(columns);
        }

        CoordinateMap(const CoordinateMap&) = delete;
        CoordinateMap& operator=(const CoordinateMap&) = delete;

        ~CoordinateMap()
        {
            delete mMap;
        }

        std::uint64_t operator()(std::uint64_t value) const
        {
            return (*mMap)(value);
        }

    private:
        ByteMap* mMap = nullptr;
    };

    std::uint32_t multiplyByX(
        std::uint32_t value,
        unsigned degree,
        std::uint64_t factor)
    {
        const bool high = (value >> (degree - 1)) & 1;
        value = (value << 1) & ((std::uint32_t{1} << degree) - 1);
        if (high)
            value ^= static_cast<std::uint32_t>(factor) & ((std::uint32_t{1} << degree) - 1);
        return value;
    }

    struct LogTable
    {
        unsigned degree = 0;
        std::uint32_t order = 0;
        std::vector<std::uint32_t> logarithm;
        std::vector<std::uint8_t> coset;
        unsigned cosetCount = 0;
    };

    LogTable makeLogTable(unsigned degree, std::uint64_t factor, std::uint32_t order)
    {
        const std::uint32_t size = std::uint32_t{1} << degree;
        LogTable table{degree, order};
        table.logarithm.assign(size, std::numeric_limits<std::uint32_t>::max());
        table.coset.assign(size, 0xff);
        for (std::uint32_t seed = 1; seed < size; ++seed)
        {
            if (table.logarithm[seed] != std::numeric_limits<std::uint32_t>::max())
                continue;
            const unsigned coset = table.cosetCount++;
            std::uint32_t current = seed;
            for (std::uint32_t exponent = 0; exponent < order; ++exponent)
            {
                if (table.logarithm[current] != std::numeric_limits<std::uint32_t>::max())
                    throw std::runtime_error("dual A4: component orbit merged early");
                table.logarithm[current] = exponent;
                table.coset[current] = static_cast<std::uint8_t>(coset);
                current = multiplyByX(current, degree, factor);
            }
            if (current != seed)
                throw std::runtime_error("dual A4: component root order mismatch");
        }
        if (table.cosetCount != ((size - 1) / order))
            throw std::runtime_error("dual A4: component coset count mismatch");
        return table;
    }

    std::uint64_t modularInverse(std::uint64_t value, std::uint64_t modulus)
    {
        if (modulus == 1)
            return 0;
        std::int64_t oldR = static_cast<std::int64_t>(value);
        std::int64_t r = static_cast<std::int64_t>(modulus);
        std::int64_t oldS = 1;
        std::int64_t s = 0;
        while (r)
        {
            const std::int64_t quotient = oldR / r;
            std::tie(oldR, r) = std::pair{r, oldR - quotient * r};
            std::tie(oldS, s) = std::pair{s, oldS - quotient * s};
        }
        if (oldR != 1)
            throw std::runtime_error("dual A4: modular inverse does not exist");
        oldS %= static_cast<std::int64_t>(modulus);
        if (oldS < 0)
            oldS += static_cast<std::int64_t>(modulus);
        return static_cast<std::uint64_t>(oldS);
    }

    struct SupportData
    {
        std::uint64_t period = 0;
        std::uint32_t orbitCount = 0;
        std::uint32_t offset = 0;
    };

    std::array<SupportData, 128> makeSupportData()
    {
        std::array<SupportData, 128> result{};
        std::uint64_t offset = 0;
        for (unsigned support = 1; support < 128; ++support)
        {
            std::uint64_t states = 1;
            std::uint64_t period = 1;
            for (unsigned component = 0; component < 7; ++component)
            {
                if ((support >> component) & 1)
                {
                    states *= (std::uint64_t{1} << Degrees[component]) - 1;
                    period = std::lcm(period, static_cast<std::uint64_t>(Orders[component]));
                }
            }
            if (states % period)
                throw std::runtime_error("dual A4: support orbit count is not integral");
            const std::uint64_t count = states / period;
            if (count > std::numeric_limits<std::uint32_t>::max() ||
                offset + count > std::numeric_limits<std::uint32_t>::max())
                throw std::runtime_error("dual A4: global orbit key overflow");
            result[support] = SupportData{
                period,
                static_cast<std::uint32_t>(count),
                static_cast<std::uint32_t>(offset)};
            offset += count;
        }
        if (offset != 204014423)
            throw std::runtime_error("dual A4: total orbit count changed");
        return result;
    }

    struct CanonicalState
    {
        std::uint32_t key = 0;
        std::uint64_t phase = 0;
        unsigned support = 0;
    };

    class Canonicalizer
    {
    public:
        Canonicalizer(
            const CoordinateMap& coordinates,
            const std::array<LogTable, 7>& logs,
            const std::array<SupportData, 128>& supports)
            : mCoordinates(coordinates), mLogs(logs), mSupports(supports)
        {
            unsigned offset = 0;
            for (unsigned component = 0; component < 7; ++component)
            {
                mOffsets[component] = offset;
                mMasks[component] =
                    ((std::uint64_t{1} << Degrees[component]) - 1) << offset;
                offset += Degrees[component];
            }
            if (offset != 64)
                throw std::runtime_error("dual A4: component coordinate width changed");
        }

        CanonicalState operator()(std::uint64_t physical) const
        {
            const std::uint64_t packed = mCoordinates(physical);
            std::array<std::uint32_t, 7> values{};
            unsigned support = 0;
            for (unsigned component = 0; component < 7; ++component)
            {
                values[component] = static_cast<std::uint32_t>(
                    (packed & mMasks[component]) >> mOffsets[component]);
                if (values[component])
                    support |= 1u << component;
            }
            if (!support)
                throw std::runtime_error("dual A4: zero pair sum encountered");

            bool initialized = false;
            std::uint64_t phase = 0;
            std::uint64_t modulus = 1;
            std::uint64_t invariant = 0;
            for (const unsigned component : CanonicalOrder)
            {
                const std::uint32_t value = values[component];
                if (!value)
                    continue;
                const LogTable& table = mLogs[component];
                const std::uint64_t log = table.logarithm[value];
                const std::uint64_t coset = table.coset[value];
                if (!initialized)
                {
                    initialized = true;
                    phase = log;
                    modulus = table.order;
                    invariant = coset;
                    continue;
                }
                invariant = invariant * table.cosetCount + coset;
                const std::uint64_t gcd = std::gcd(modulus, static_cast<std::uint64_t>(table.order));
                const std::uint64_t phaseModGcd = phase % gcd;
                const std::uint64_t discrepancy = (log + gcd - phaseModGcd) % gcd;
                invariant = invariant * gcd + discrepancy;
                const std::uint64_t adjusted = (log + table.order - discrepancy) % table.order;
                const std::uint64_t reducedRight = table.order / gcd;
                std::uint64_t multiple = 0;
                if (reducedRight != 1)
                {
                    const std::int64_t difference =
                        static_cast<std::int64_t>(adjusted) -
                        static_cast<std::int64_t>(phase % table.order);
                    std::int64_t quotient = difference / static_cast<std::int64_t>(gcd);
                    quotient %= static_cast<std::int64_t>(reducedRight);
                    if (quotient < 0)
                        quotient += static_cast<std::int64_t>(reducedRight);
                    const std::uint64_t inverse = modularInverse(
                        (modulus / gcd) % reducedRight, reducedRight);
                    multiple = (static_cast<std::uint64_t>(quotient) * inverse) % reducedRight;
                }
                phase += modulus * multiple;
                modulus *= reducedRight;
                phase %= modulus;
            }
            const SupportData& data = mSupports[support];
            if (modulus != data.period || invariant >= data.orbitCount)
                throw std::runtime_error("dual A4: canonical quotient mismatch");
            return CanonicalState{
                static_cast<std::uint32_t>(data.offset + invariant), phase, support};
        }

    private:
        const CoordinateMap& mCoordinates;
        const std::array<LogTable, 7>& mLogs;
        const std::array<SupportData, 128>& mSupports;
        std::array<unsigned, 7> mOffsets{};
        std::array<std::uint64_t, 7> mMasks{};
    };

    struct Record
    {
        std::uint64_t phase = 0;
        std::uint32_t key = 0;
        std::uint16_t length = 0;
        std::uint16_t reserved = 0;
    };
    static_assert(sizeof(Record) == 16);

    struct Event
    {
        std::uint64_t position = 0;
        std::int32_t delta = 0;
    };

    U128 chooseTwo(std::uint64_t value)
    {
        return (value & 1)
            ? multiply(value, (value - 1) / 2)
            : multiply(value / 2, value - 1);
    }

    void addMultiplicity(
        std::vector<U128>& histogram,
        std::uint64_t multiplicity,
        std::uint64_t stateCount)
    {
        if (!multiplicity || !stateCount)
            return;
        if (multiplicity >= histogram.size())
            histogram.resize(static_cast<std::size_t>(multiplicity) + 1);
        histogram[multiplicity] += U128{stateCount};
    }

    U128 collisionCount(
        const std::vector<Record>& records,
        std::size_t first,
        std::size_t last,
        std::uint64_t period,
        std::vector<Event>& events,
        std::vector<U128>& multiplicityHistogram)
    {
#if INDEPENDENT_VERIFIER
        struct Arc
        {
            std::array<std::pair<std::uint64_t, std::uint64_t>, 2> segments{};
            unsigned count = 0;
        };
        auto arc = [period](const Record& record)
        {
            Arc result;
            const std::uint64_t remainder = record.length % period;
            if (!remainder)
                return result;
            const std::uint64_t end = record.phase + remainder;
            if (end <= period)
                result.segments[result.count++] = {record.phase, end};
            else
            {
                result.segments[result.count++] = {record.phase, period};
                result.segments[result.count++] = {0, end - period};
            }
            return result;
        };
        auto overlap = [](const Arc& left, const Arc& right)
        {
            std::uint64_t result = 0;
            for (unsigned i = 0; i < left.count; ++i)
            {
                for (unsigned j = 0; j < right.count; ++j)
                {
                    const std::uint64_t begin = std::max(
                        left.segments[i].first, right.segments[j].first);
                    const std::uint64_t end = std::min(
                        left.segments[i].second, right.segments[j].second);
                    if (begin < end)
                        result += end - begin;
                }
            }
            return result;
        };
        U128 result;
        for (std::size_t left = first; left < last; ++left)
        {
            const std::uint64_t leftQuotient = records[left].length / period;
            const std::uint64_t leftRemainder = records[left].length % period;
            result += chooseTwo(leftQuotient) * period;
            result += U128{leftQuotient} * leftRemainder;
            const Arc leftArc = arc(records[left]);
            for (std::size_t right = left + 1; right < last; ++right)
            {
                const std::uint64_t rightQuotient = records[right].length / period;
                const std::uint64_t rightRemainder = records[right].length % period;
                result += multiply(leftQuotient, rightQuotient) * period;
                result += U128{leftQuotient} * rightRemainder;
                result += U128{rightQuotient} * leftRemainder;
                result += U128{overlap(leftArc, arc(records[right]))};
            }
        }
        std::vector<std::uint64_t> breakpoints{0, period};
        breakpoints.reserve(4 * (last - first) + 2);
        for (std::size_t index = first; index < last; ++index)
        {
            const Arc current = arc(records[index]);
            for (unsigned segment = 0; segment < current.count; ++segment)
            {
                breakpoints.push_back(current.segments[segment].first);
                breakpoints.push_back(current.segments[segment].second);
            }
        }
        std::sort(breakpoints.begin(), breakpoints.end());
        breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()), breakpoints.end());
        U128 histogramCollision;
        for (std::size_t point = 0; point + 1 < breakpoints.size(); ++point)
        {
            const std::uint64_t begin = breakpoints[point];
            const std::uint64_t length = breakpoints[point + 1] - begin;
            if (!length)
                continue;
            std::uint64_t multiplicity = 0;
            for (std::size_t index = first; index < last; ++index)
            {
                multiplicity += records[index].length / period;
                const Arc current = arc(records[index]);
                for (unsigned segment = 0; segment < current.count; ++segment)
                {
                    if (current.segments[segment].first <= begin &&
                        begin < current.segments[segment].second)
                    {
                        ++multiplicity;
                        break;
                    }
                }
            }
            addMultiplicity(multiplicityHistogram, multiplicity, length);
            histogramCollision += chooseTwo(multiplicity) * length;
        }
        if (!(histogramCollision == result))
            throw std::runtime_error("dual A4: independent multiplicity histogram differs");
        return result;
#else
        std::uint64_t uniform = 0;
        std::uint64_t remainderMass = 0;
        events.clear();
        events.reserve(4 * (last - first));
        auto addArc = [&](std::uint64_t begin, std::uint64_t end)
        {
            if (begin == end)
                return;
            events.push_back(Event{begin, 1});
            events.push_back(Event{end, -1});
        };
        for (std::size_t index = first; index < last; ++index)
        {
            const Record& record = records[index];
            uniform += record.length / period;
            const std::uint64_t remainder = record.length % period;
            remainderMass += remainder;
            if (!remainder)
                continue;
            const std::uint64_t end = record.phase + remainder;
            if (end <= period)
                addArc(record.phase, end);
            else
            {
                addArc(record.phase, period);
                addArc(0, end - period);
            }
        }
        U128 result = chooseTwo(uniform) * period + U128{uniform} * remainderMass;
        if (events.empty())
        {
            addMultiplicity(multiplicityHistogram, uniform, period);
            return result;
        }
        std::sort(events.begin(), events.end(), [](const Event& left, const Event& right)
        {
            return left.position < right.position;
        });
        std::int64_t active = 0;
        std::uint64_t previous = 0;
        std::size_t index = 0;
        while (index < events.size())
        {
            const std::uint64_t position = events[index].position;
            if (active < 0)
                throw std::runtime_error("dual A4: negative interval coverage");
            addMultiplicity(
                multiplicityHistogram,
                uniform + static_cast<std::uint64_t>(active),
                position - previous);
            result += chooseTwo(static_cast<std::uint64_t>(active)) * (position - previous);
            std::int64_t delta = 0;
            while (index < events.size() && events[index].position == position)
                delta += events[index++].delta;
            active += delta;
            previous = position;
        }
        if (active != 0 || previous > period)
            throw std::runtime_error("dual A4: cyclic interval sweep did not close");
        addMultiplicity(
            multiplicityHistogram,
            uniform + static_cast<std::uint64_t>(active),
            period - previous);
        result += chooseTwo(static_cast<std::uint64_t>(active)) * (period - previous);
        return result;
#endif
    }
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("usage: dual_a4.exe OUTPUT.json");
        const auto started = std::chrono::steady_clock::now();
        const auto parityColumns = systematicRightColumns();
        const ByteMap parity(parityColumns);
        std::array<std::uint64_t, 64> outputStepColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            outputStepColumns[bit] = accumulate(parity(std::uint64_t{1} << bit));
        const ByteMap transposeStep(transposeColumns(outputStepColumns));

        std::array<std::vector<std::uint64_t>, 7> componentKernels;
        std::vector<std::uint64_t> krylovBasis;
        for (unsigned component = 0; component < 7; ++component)
        {
            componentKernels[component] = kernelBasis(
                transposeStep, Factors[component], Degrees[component]);
#if INDEPENDENT_VERIFIER
            std::uint64_t current = componentKernels[component].back();
#else
            std::uint64_t current = componentKernels[component][0];
#endif
            for (unsigned exponent = 0; exponent < Degrees[component]; ++exponent)
            {
                krylovBasis.push_back(current);
                current = transposeStep(current);
            }
        }
        const CoordinateMap coordinateMap(krylovBasis);
        std::array<LogTable, 7> logTables;
        for (unsigned component = 0; component < 7; ++component)
            logTables[component] = makeLogTable(
                Degrees[component], Factors[component], Orders[component]);
        const auto supports = makeSupportData();
        const Canonicalizer canonicalize(coordinateMap, logTables, supports);

        for (unsigned bit = 0; bit < 64; ++bit)
        {
            const CanonicalState left = canonicalize(std::uint64_t{1} << bit);
            const CanonicalState right = canonicalize(transposeStep(std::uint64_t{1} << bit));
            if (left.key != right.key ||
                right.phase != (left.phase + 1) % supports[left.support].period)
                throw std::runtime_error("dual A4: canonical phase self-check failed");
        }

        std::vector<Record> records;
        records.reserve(ExpectedRecordCount);
        std::array<std::uint64_t, 64> shifted{};
        for (unsigned bit = 0; bit < 64; ++bit)
            shifted[bit] = std::uint64_t{1} << bit;
        auto append = [&](std::uint64_t state, unsigned length)
        {
            const CanonicalState canonical = canonicalize(state);
            records.push_back(Record{
                canonical.phase,
                canonical.key,
                static_cast<std::uint16_t>(length),
                0});
        };
        for (unsigned left = 0; left < 64; ++left)
        {
            for (unsigned right = left + 1; right < 64; ++right)
                append((std::uint64_t{1} << left) ^ (std::uint64_t{1} << right), Nodes);
        }
        for (unsigned gap = 1; gap < Nodes; ++gap)
        {
            for (unsigned right = 0; right < 64; ++right)
                shifted[right] = transposeStep(shifted[right]);
            const unsigned length = Nodes - gap;
            for (unsigned left = 0; left < 64; ++left)
            {
                const std::uint64_t first = std::uint64_t{1} << left;
                for (unsigned right = 0; right < 64; ++right)
                    append(first ^ shifted[right], length);
            }
            if ((gap & 4095) == 0)
                std::cerr << "generated_gap=" << gap << " records=" << records.size() << '\n';
        }
        if (records.size() != ExpectedRecordCount)
            throw std::runtime_error("dual A4: normalized pair record count mismatch");
        const double generationSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::cerr << "generation_complete records=" << records.size()
                  << " seconds=" << generationSeconds << '\n';

        const auto sortingStarted = std::chrono::steady_clock::now();
        std::sort(records.begin(), records.end(), [](const Record& left, const Record& right)
        {
            if (left.key != right.key)
                return left.key < right.key;
            if (left.phase != right.phase)
                return left.phase < right.phase;
            return left.length < right.length;
        });
        const double sortingSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - sortingStarted).count();
        std::cerr << "sorting_complete seconds=" << sortingSeconds << '\n';

        U128 pairCollisions = 0;
        std::uint64_t occupiedOrbitKeys = 0;
        std::uint64_t maximumRecordsPerOrbit = 0;
        std::vector<Event> events;
        std::vector<U128> multiplicityHistogram;
        unsigned support = 1;
        for (std::size_t first = 0; first < records.size();)
        {
            std::size_t last = first + 1;
            while (last < records.size() && records[last].key == records[first].key)
                ++last;
            while (support < 128 &&
                records[first].key >=
                    std::uint64_t{supports[support].offset} + supports[support].orbitCount)
                ++support;
            if (support == 128 || records[first].key < supports[support].offset)
                throw std::runtime_error("dual A4: orbit key has no support range");
            pairCollisions += collisionCount(
                records,
                first,
                last,
                supports[support].period,
                events,
                multiplicityHistogram);
            ++occupiedOrbitKeys;
            maximumRecordsPerOrbit = std::max<std::uint64_t>(
                maximumRecordsPerOrbit, last - first);
            first = last;
        }
        if (remainder(pairCollisions, 3))
            throw std::runtime_error("dual A4: pair collisions are not divisible by three");
        unsigned divisionRemainder = 0;
        const U128 dualA4 = divide(pairCollisions, 3, divisionRemainder);
        if (divisionRemainder)
            throw std::runtime_error("dual A4: division by three left a remainder");
        U128 distinctPairSums;
        U128 pairMass;
        U128 histogramCollisions;
        std::uint64_t maximumPairSumMultiplicity = 0;
        for (std::size_t multiplicity = 1; multiplicity < multiplicityHistogram.size(); ++multiplicity)
        {
            if (multiplicityHistogram[multiplicity].low == 0 &&
                multiplicityHistogram[multiplicity].high == 0)
                continue;
            maximumPairSumMultiplicity = multiplicity;
            distinctPairSums += multiplicityHistogram[multiplicity];
            pairMass += multiplicityHistogram[multiplicity] * multiplicity;
            histogramCollisions += multiplicityHistogram[multiplicity] *
                chooseTwo(multiplicity).low;
        }
        const std::uint64_t responseCoordinates = std::uint64_t{Nodes} * Coordinates;
        const U128 expectedPairMass = multiply(
            responseCoordinates, responseCoordinates - 1) * 1;
        unsigned expectedPairMassRemainder = 0;
        const U128 expectedUnorderedPairs = divide(
            expectedPairMass, 2, expectedPairMassRemainder);
        if (expectedPairMassRemainder || !(pairMass == expectedUnorderedPairs))
            throw std::runtime_error("dual A4: pair multiplicity mass mismatch");
        if (!(histogramCollisions == pairCollisions))
            throw std::runtime_error("dual A4: pair multiplicity collision mismatch");
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        std::ofstream output(argv[1]);
        if (!output)
            throw std::runtime_error("dual A4: cannot open output receipt");
        output << "{\n"
               << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal02-dual-a4-v1\",\n"
               << "  \"candidate\": \"Riffle PacketMul-WrapMul-2Lap g=4\",\n"
               << "  \"evidence_label\": \""
               << (INDEPENDENT_VERIFIER
                    ? "INDEPENDENT_EXACT_PAIRWISE_ARC_REPLAY"
                    : "EXACT_NORMALIZED_PAIR_ORBIT_CENSUS")
               << "\",\n"
               << "  \"response_nodes\": " << Nodes << ",\n"
               << "  \"response_coordinates\": " << std::uint64_t{Nodes} * Coordinates << ",\n"
               << "  \"irreducible_component_degrees\": [1, 2, 4, 9, 10, 18, 20],\n"
               << "  \"irreducible_component_orders\": [1, 3, 15, 511, 1023, 37449, 1048575],\n"
               << "  \"total_autonomous_orbit_keys\": 204014423,\n"
               << "  \"normalized_unordered_pair_types\": " << records.size() << ",\n"
               << "  \"occupied_pair_sum_orbit_keys\": " << occupiedOrbitKeys << ",\n"
               << "  \"maximum_normalized_pair_types_per_orbit_key\": " << maximumRecordsPerOrbit << ",\n"
               << "  \"distinct_nonzero_pair_sums\": \"" << decimal(distinctPairSums) << "\",\n"
               << "  \"maximum_pair_sum_multiplicity\": " << maximumPairSumMultiplicity << ",\n"
               << "  \"pair_sum_multiplicity_histogram\": [";
        bool firstHistogram = true;
        for (std::size_t multiplicity = 1; multiplicity < multiplicityHistogram.size(); ++multiplicity)
        {
            if (multiplicityHistogram[multiplicity].low == 0 &&
                multiplicityHistogram[multiplicity].high == 0)
                continue;
            output << (firstHistogram ? "\n" : ",\n")
                   << "    {\"multiplicity\": " << multiplicity
                   << ", \"state_count\": \"" << decimal(multiplicityHistogram[multiplicity]) << "\"}";
            firstHistogram = false;
        }
        if (!firstHistogram)
            output << '\n';
        output << "  ],\n"
               << "  \"equal_pair_sum_collisions\": \"" << decimal(pairCollisions) << "\",\n"
               << "  \"dual_weight4_word_count\": \"" << decimal(dualA4) << "\",\n"
               << "  \"collision_to_word_divisor\": 3,\n"
               << "  \"generation_seconds\": " << generationSeconds << ",\n"
               << "  \"sorting_seconds\": " << sortingSeconds << ",\n"
               << "  \"elapsed_seconds\": " << elapsed << ",\n"
               << "  \"result\": \"PASS\"\n"
               << "}\n";
        if (!output)
            throw std::runtime_error("dual A4: receipt write failed");
        std::cout << "normalized_pair_types=" << records.size() << '\n'
                  << "occupied_orbit_keys=" << occupiedOrbitKeys << '\n'
                  << "maximum_pair_sum_multiplicity=" << maximumPairSumMultiplicity << '\n'
                  << "pair_collisions=" << decimal(pairCollisions) << '\n'
                  << "dual_A4=" << decimal(dualA4) << '\n'
                  << "elapsed_seconds=" << elapsed << '\n'
                  << "output=" << argv[1] << '\n'
                  << "status=EXACT_NORMALIZED_PAIR_ORBIT_CENSUS\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
