# Riffle FieldCheckpointAccumulate t=32 s=256 K=32

This exploratory candidate uses 256 accumulator lanes and retains a 1024-bit
checkpoint epoch. Each lane is visited four times per epoch.

Under the modeled outer spectrum, the exact floating-point calculation for
512 through 4096 regular active blocks has 21609.2409 bits of aggregate margin
at 9% relative distance. All-one outer words and outward rounding remain open.

The transposed mixing cost is one width-256 field map per 1024 output bits.
That cost has not been benchmarked.
