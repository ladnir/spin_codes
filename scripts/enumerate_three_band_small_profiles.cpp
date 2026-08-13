#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Partition {
    std::array<std::uint8_t, 6> label{};
    std::array<std::uint8_t, 6> clusterSize{};
    std::uint16_t edgeMask = 0;
};

void generatePartitions(int size, int index, int maximum,
                        std::array<std::uint8_t, 6>& labels,
                        std::vector<Partition>& result) {
    if (index == size) {
        Partition partition;
        partition.label = labels;
        std::array<std::uint8_t, 6> counts{};
        for (int vertex = 0; vertex < size; ++vertex) {
            ++counts[labels[vertex]];
        }
        for (int vertex = 0; vertex < size; ++vertex) {
            partition.clusterSize[vertex] = counts[labels[vertex]];
        }
        int bit = 0;
        for (int left = 0; left < size; ++left) {
            for (int right = left + 1; right < size; ++right, ++bit) {
                if (labels[left] == labels[right]) {
                    partition.edgeMask |= static_cast<std::uint16_t>(1u << bit);
                }
            }
        }
        result.push_back(partition);
        return;
    }
    for (int label = 0; label <= maximum + 1; ++label) {
        labels[index] = static_cast<std::uint8_t>(label);
        generatePartitions(size, index + 1, std::max(maximum, label), labels,
                           result);
    }
}

int componentCount(int size, std::uint16_t mask) {
    std::array<int, 6> parent{};
    for (int vertex = 0; vertex < size; ++vertex) parent[vertex] = vertex;
    const auto find = [&](int vertex, auto&& self) -> int {
        if (parent[vertex] != vertex) parent[vertex] = self(parent[vertex], self);
        return parent[vertex];
    };
    int bit = 0;
    for (int left = 0; left < size; ++left) {
        for (int right = left + 1; right < size; ++right, ++bit) {
            if ((mask & static_cast<std::uint16_t>(1u << bit)) == 0) continue;
            const int leftRoot = find(left, find);
            const int rightRoot = find(right, find);
            if (leftRoot != rightRoot) parent[rightRoot] = leftRoot;
        }
    }
    int components = 0;
    for (int vertex = 0; vertex < size; ++vertex) {
        if (find(vertex, find) == vertex) ++components;
    }
    return components;
}

std::uint64_t profileKey(int size, int components,
                         const Partition& first, const Partition& second,
                         const Partition& third) {
    std::array<std::uint16_t, 6> kinds{};
    for (int vertex = 0; vertex < size; ++vertex) {
        kinds[vertex] = static_cast<std::uint16_t>(
            first.clusterSize[vertex]
            | (second.clusterSize[vertex] << 3)
            | (third.clusterSize[vertex] << 6));
    }
    std::sort(kinds.begin(), kinds.begin() + size);
    std::uint64_t key = static_cast<std::uint64_t>(components);
    int shift = 3;
    for (int vertex = 0; vertex < size; ++vertex, shift += 9) {
        key |= static_cast<std::uint64_t>(kinds[vertex]) << shift;
    }
    return key;
}

void enumerateSize(int size) {
    std::array<std::uint8_t, 6> labels{};
    std::vector<Partition> partitions;
    labels[0] = 0;
    generatePartitions(size, 1, 0, labels, partitions);

    std::array<std::uint8_t, 1u << 15> components{};
    const int edgeBits = size * (size - 1) / 2;
    for (int mask = 0; mask < (1 << edgeBits); ++mask) {
        components[mask] = static_cast<std::uint8_t>(
            componentCount(size, static_cast<std::uint16_t>(mask)));
    }

    std::unordered_map<std::uint64_t, std::uint64_t> profiles;
    std::uint64_t compatibleTriples = 0;
    for (const Partition& first : partitions) {
        for (const Partition& second : partitions) {
            if (first.edgeMask & second.edgeMask) continue;
            for (const Partition& third : partitions) {
                if ((first.edgeMask & third.edgeMask) ||
                    (second.edgeMask & third.edgeMask)) {
                    continue;
                }
                const auto joined = static_cast<std::uint16_t>(
                    first.edgeMask | second.edgeMask | third.edgeMask);
                const std::uint64_t key = profileKey(
                    size, components[joined], first, second, third);
                ++profiles[key];
                ++compatibleTriples;
            }
        }
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> ordered(
        profiles.begin(), profiles.end());
    std::sort(ordered.begin(), ordered.end());
    std::cout << "summary " << size << ' ' << partitions.size() << ' '
              << compatibleTriples << ' ' << profiles.size() << '\n';
    for (const auto& [key, count] : ordered) {
        const int componentsValue = static_cast<int>(key & 7u);
        std::cout << "profile " << size << ' ' << componentsValue << ' '
                  << count;
        int shift = 3;
        for (int vertex = 0; vertex < size; ++vertex, shift += 9) {
            const int kind = static_cast<int>((key >> shift) & 0x1ffu);
            const int firstDegree = kind & 7;
            const int secondDegree = (kind >> 3) & 7;
            const int thirdDegree = (kind >> 6) & 7;
            std::cout << ' ' << firstDegree << ',' << secondDegree << ','
                      << thirdDegree;
        }
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    for (int size = 2; size <= 6; ++size) enumerateSize(size);
    std::cout << "status EXACT_PACKED_PROFILE_ENUMERATION\n";
    return 0;
}
