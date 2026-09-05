# Riffle FieldCheckpointAccumulate t=32 s=128 K=16

This exploratory candidate uses 128 accumulator lanes and refreshes the state
every 512 output bits. Each lane is visited four times per epoch.

Under the modeled outer spectrum, the exact floating-point calculation for
512 through 4096 regular active blocks has 21135.6506 bits of aggregate margin
at 9% relative distance. All-one outer words and outward rounding remain open.

The transposed mixing cost is two width-128 field maps per 1024 output bits.
That cost has not been benchmarked.
