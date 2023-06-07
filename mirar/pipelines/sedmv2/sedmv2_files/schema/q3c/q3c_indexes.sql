CREATE EXTENSION q3c;

CREATE INDEX IF NOT EXISTS raw_q3c_ang2ipix_idx ON raw (q3c_ang2ipix(ra, dec));
CLUSTER raw_q3c_ang2ipix_idx ON raw;
ANALYZE raw;

CREATE INDEX IF NOT EXISTS proc_q3c_ang2ipix_idx ON proc (q3c_ang2ipix(crval1, crval2));
CLUSTER proc_q3c_ang2ipix_idx ON proc;
ANALYZE proc;
