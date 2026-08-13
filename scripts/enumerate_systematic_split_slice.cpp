// Exact C[t,q] rows for the systematic-left split of the committed
// extended-BCH [128,64,22] code.  Here t is the left-half/input weight and q
// is the right-half/state weight.  One enumeration also produces 64-t.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

namespace
{
	constexpr std::uint64_t generator = 0xf4845518b9582a1fULL;
	std::array<std::uint64_t, 64> stateRows{};
	std::uint64_t allState = 0;
	using Histogram = std::array<std::uint64_t, 65>;

	struct Task
	{
		unsigned start;
		unsigned remaining;
		std::uint64_t state;
	};

	std::uint64_t applyLinear(
		const std::array<std::uint64_t, 64>& columns,
		std::uint64_t value)
	{
		std::uint64_t result = 0;
		while (value)
		{
			const auto bit = value & (0 - value);
			result ^= columns[std::countr_zero(bit)];
			value ^= bit;
		}
		return result;
	}

	std::array<std::uint64_t, 64> invertLinear(
		const std::array<std::uint64_t, 64>& columns)
	{
		std::array<std::uint64_t, 64> matrixRows{};
		std::array<std::uint64_t, 64> inverseRows{};
		for (unsigned row = 0; row < 64; ++row)
		{
			for (unsigned column = 0; column < 64; ++column)
				matrixRows[row] |= ((columns[column] >> row) & 1) << column;
			inverseRows[row] = std::uint64_t{ 1 } << row;
		}
		for (unsigned column = 0; column < 64; ++column)
		{
			const auto pivot = std::find_if(
				matrixRows.begin() + column,
				matrixRows.end(),
				[column](std::uint64_t row) { return (row >> column) & 1; });
			if (pivot == matrixRows.end())
				std::abort();
			const auto pivotIndex = static_cast<unsigned>(pivot - matrixRows.begin());
			std::swap(matrixRows[column], matrixRows[pivotIndex]);
			std::swap(inverseRows[column], inverseRows[pivotIndex]);
			for (unsigned row = 0; row < 64; ++row)
				if (row != column && ((matrixRows[row] >> column) & 1))
				{
					matrixRows[row] ^= matrixRows[column];
					inverseRows[row] ^= inverseRows[column];
				}
		}
		std::array<std::uint64_t, 64> inverseColumns{};
		for (unsigned output = 0; output < 64; ++output)
			for (unsigned input = 0; input < 64; ++input)
				inverseColumns[output] |= ((inverseRows[input] >> output) & 1) << input;
		return inverseColumns;
	}

	void initializeSystematicRows()
	{
		std::array<std::uint64_t, 64> leftColumns{};
		std::array<std::uint64_t, 64> rightColumns{};
		for (unsigned column = 0; column < 64; ++column)
		{
			leftColumns[column] = generator << column;
			rightColumns[column] =
				(column ? generator >> (64 - column) : 0) |
				(std::uint64_t{ 1 } << 63);
		}

		std::array<std::uint64_t, 64> matrixRows{};
		std::array<std::uint64_t, 64> inverseRows{};
		for (unsigned row = 0; row < 64; ++row)
		{
			for (unsigned column = 0; column < 64; ++column)
				matrixRows[row] |= ((leftColumns[column] >> row) & 1) << column;
			inverseRows[row] = std::uint64_t{ 1 } << row;
		}
		for (unsigned column = 0; column < 64; ++column)
		{
			const auto pivot = std::find_if(
				matrixRows.begin() + column,
				matrixRows.end(),
				[column](std::uint64_t row) { return (row >> column) & 1; });
			if (pivot == matrixRows.end())
				std::abort();
			const auto pivotIndex = static_cast<unsigned>(pivot - matrixRows.begin());
			std::swap(matrixRows[column], matrixRows[pivotIndex]);
			std::swap(inverseRows[column], inverseRows[pivotIndex]);
			for (unsigned row = 0; row < 64; ++row)
				if (row != column && ((matrixRows[row] >> column) & 1))
				{
					matrixRows[row] ^= matrixRows[column];
					inverseRows[row] ^= inverseRows[column];
				}
		}

		std::array<std::uint64_t, 64> inverseColumns{};
		for (unsigned output = 0; output < 64; ++output)
			for (unsigned input = 0; input < 64; ++input)
				inverseColumns[output] |= ((inverseRows[input] >> output) & 1) << input;
		for (unsigned column = 0; column < 64; ++column)
		{
			if (applyLinear(leftColumns, inverseColumns[column]) != (std::uint64_t{ 1 } << column))
				std::abort();
			stateRows[column] = applyLinear(rightColumns, inverseColumns[column]);
			allState ^= stateRows[column];
		}
		if (allState != ~std::uint64_t{ 0 })
			std::abort();
	}

	void enumerate(
		unsigned start,
		unsigned remaining,
		std::uint64_t state,
		Histogram& lowHistogram,
		Histogram& highHistogram)
	{
		if (remaining == 0)
		{
			++lowHistogram[std::popcount(state)];
			++highHistogram[std::popcount(state ^ allState)];
			return;
		}
		const auto last = 64u - remaining;
		for (auto index = start; index <= last; ++index)
			enumerate(
				index + 1,
				remaining - 1,
				state ^ stateRows[index],
				lowHistogram,
				highHistogram);
	}

	void print(unsigned inputWeight, const Histogram& histogram)
	{
		for (unsigned stateWeight = 0; stateWeight <= 64; ++stateWeight)
			if (histogram[stateWeight])
				std::cout << inputWeight << ',' << stateWeight << ','
					<< histogram[stateWeight] << '\n';
	}

	void printColumn(unsigned stateWeight, const Histogram& histogram)
	{
		for (unsigned inputWeight = 0; inputWeight <= 64; ++inputWeight)
			if (histogram[inputWeight])
				std::cout << inputWeight << ',' << stateWeight << ','
					<< histogram[inputWeight] << '\n';
	}
}

int main(int argc, char** argv)
{
	const auto weight = argc > 1 ? static_cast<unsigned>(std::strtoul(argv[1], nullptr, 0)) : 7;
	const auto requestedThreads = argc > 2
		? static_cast<unsigned>(std::strtoul(argv[2], nullptr, 0))
		: std::thread::hardware_concurrency();
	const auto columnMode = argc > 3 && std::string_view(argv[3]) == "column";
	if (weight > 32)
	{
		std::cerr << "weight must be at most 32\n";
		return 1;
	}
	initializeSystematicRows();
	if (columnMode)
	{
		stateRows = invertLinear(stateRows);
		allState = 0;
		for (const auto row : stateRows)
			allState ^= row;
		if (allState != ~std::uint64_t{ 0 })
			std::abort();
	}

	std::vector<Task> tasks;
	if (weight == 0)
		tasks.push_back({ 0, 0, 0 });
	else if (weight == 1)
		for (unsigned first = 0; first < 64; ++first)
			tasks.push_back({ first + 1, 0, stateRows[first] });
	else
		for (unsigned first = 0; first + weight <= 64; ++first)
			for (unsigned second = first + 1; second + weight - 1 <= 64; ++second)
				tasks.push_back({
					second + 1,
					weight - 2,
					stateRows[first] ^ stateRows[second]
				});

	const auto threadCount = std::max(
		1u,
		std::min(
			requestedThreads ? requestedThreads : 1u,
			static_cast<unsigned>(tasks.size())));
	std::vector<Histogram> lowHistograms(threadCount);
	std::vector<Histogram> highHistograms(threadCount);
	std::atomic<std::size_t> nextTask = 0;
	std::vector<std::thread> workers;
	const auto started = std::chrono::steady_clock::now();
	for (unsigned threadIndex = 0; threadIndex < threadCount; ++threadIndex)
		workers.emplace_back([&, threadIndex] {
			while (true)
			{
				const auto taskIndex = nextTask.fetch_add(1, std::memory_order_relaxed);
				if (taskIndex >= tasks.size())
					break;
				const auto& task = tasks[taskIndex];
				enumerate(
					task.start,
					task.remaining,
					task.state,
					lowHistograms[threadIndex],
					highHistograms[threadIndex]);
			}
		});
	for (auto& worker : workers)
		worker.join();

	Histogram lowHistogram{};
	Histogram highHistogram{};
	for (unsigned threadIndex = 0; threadIndex < threadCount; ++threadIndex)
		for (unsigned stateWeight = 0; stateWeight <= 64; ++stateWeight)
		{
			lowHistogram[stateWeight] += lowHistograms[threadIndex][stateWeight];
			highHistogram[stateWeight] += highHistograms[threadIndex][stateWeight];
		}
	const auto elapsed = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - started).count();
	std::cerr << "enumerated C(64," << weight << ") with " << threadCount
		<< " threads in " << elapsed << " seconds\n";
	std::cout << "input_weight,state_weight,count\n";
	if (columnMode)
		printColumn(weight, lowHistogram);
	else
		print(weight, lowHistogram);
	if (weight != 32)
	{
		if (columnMode)
			printColumn(64 - weight, highHistogram);
		else
			print(64 - weight, highHistogram);
	}
}
