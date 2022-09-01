CREATE EXTENSION q3c;

CREATE INDEX IF NOT EXISTS ON raw (q3c_ang2ipix(ra, dec));
CLUSTER raw_q3c_ang2ipix_idx ON raw;
ANALYZE raw;