PRAGMA foreign_keys = ON;

CREATE TABLE schema_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE sources (
    source_id INTEGER PRIMARY KEY,
    path TEXT NOT NULL UNIQUE,
    sha256 TEXT NOT NULL,
    ingestor TEXT NOT NULL,
    source_status TEXT NOT NULL,
    transfer_review_status TEXT NOT NULL DEFAULT 'historical_pending_review',
    notes TEXT NOT NULL DEFAULT ''
);

CREATE TABLE outer_models (
    outer_id TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    family TEXT NOT NULL CHECK (family IN ('bch', 'rm', 'random')),
    model_kind TEXT NOT NULL CHECK (
        model_kind IN (
            'fixed_exact_spectrum',
            'fixed_certified_constraints',
            'fixed_estimated_spectrum',
            'random_ensemble_expectation'
        )
    ),
    block_bits INTEGER NOT NULL,
    dimension INTEGER NOT NULL,
    minimum_distance INTEGER,
    reuse_policy TEXT NOT NULL,
    UNIQUE(label, model_kind)
);

CREATE TABLE inner_configs (
    inner_id TEXT PRIMARY KEY,
    family TEXT NOT NULL CHECK (family = 'rm2sub'),
    step_bits INTEGER NOT NULL,
    state_bits INTEGER NOT NULL,
    persistence_exponent REAL NOT NULL,
    map_tag TEXT NOT NULL,
    a_minimum_distance INTEGER,
    kernel_minimum_distance INTEGER,
    kernel_weight_four INTEGER,
    UNIQUE(step_bits, state_bits, map_tag)
);

CREATE TABLE results (
    result_id TEXT PRIMARY KEY,
    source_id INTEGER NOT NULL REFERENCES sources(source_id),
    source_locator TEXT NOT NULL,
    study TEXT NOT NULL,
    outer_id TEXT NOT NULL REFERENCES outer_models(outer_id),
    inner_id TEXT NOT NULL REFERENCES inner_configs(inner_id),
    message_exponent INTEGER,
    message_bits INTEGER NOT NULL,
    output_bits INTEGER NOT NULL,
    outer_rows INTEGER NOT NULL,
    epochs_per_region INTEGER,
    bad_weight INTEGER NOT NULL,
    distance_target REAL NOT NULL,
    occupation_min INTEGER NOT NULL,
    occupation_max INTEGER NOT NULL,
    coverage_kind TEXT NOT NULL CHECK (
        coverage_kind IN ('single_occupation', 'occupation_union', 'full_distance')
    ),
    result_class TEXT NOT NULL CHECK (
        result_class IN ('certified', 'diagnostic', 'estimated', 'reference')
    ),
    arithmetic TEXT NOT NULL,
    margin_bits REAL NOT NULL,
    margin_bits_text TEXT,
    failure_upper_text TEXT,
    dominant_weight INTEGER,
    dominant_witness_at_grid_edge INTEGER,
    comparison_eligible INTEGER NOT NULL DEFAULT 0 CHECK (comparison_eligible IN (0,1)),
    notes TEXT NOT NULL DEFAULT '',
    UNIQUE(source_id, source_locator)
);

CREATE INDEX results_parameter_index
ON results(message_bits, outer_id, inner_id, distance_target);

CREATE INDEX results_occupation_index
ON results(occupation_min, occupation_max, result_class);

CREATE VIEW landscape AS
SELECT
    r.result_id,
    r.study,
    r.outer_id,
    r.inner_id,
    o.family AS outer_family,
    o.label AS outer_label,
    o.model_kind AS outer_model_kind,
    o.block_bits,
    o.dimension,
    o.minimum_distance,
    i.step_bits,
    i.state_bits,
    i.persistence_exponent,
    i.persistence_exponent - COALESCE(r.message_exponent, 0) AS persistence_offset,
    i.map_tag,
    r.message_exponent,
    r.message_bits,
    r.output_bits,
    r.outer_rows,
    r.epochs_per_region,
    r.bad_weight,
    r.distance_target,
    r.occupation_min,
    r.occupation_max,
    r.coverage_kind,
    r.result_class,
    r.arithmetic,
    r.margin_bits,
    r.margin_bits_text,
    r.failure_upper_text,
    r.dominant_weight,
    r.dominant_witness_at_grid_edge,
    r.comparison_eligible,
    s.path AS source_path,
    s.sha256 AS source_sha256,
    s.transfer_review_status,
    r.source_locator,
    r.notes
FROM results AS r
JOIN outer_models AS o USING (outer_id)
JOIN inner_configs AS i USING (inner_id)
JOIN sources AS s USING (source_id);

CREATE VIEW q1_curves AS
SELECT * FROM landscape
WHERE occupation_min = 1
  AND occupation_max = 1
  AND coverage_kind = 'single_occupation';

CREATE VIEW certified_results AS
SELECT * FROM landscape
WHERE result_class = 'certified' AND transfer_review_status = 'activation_aware';

CREATE VIEW results_under_review AS
SELECT * FROM landscape WHERE transfer_review_status <> 'activation_aware';
