CREATE EXTENSION q3c;

CREATE INDEX ON raw (q3c_ang2pix(ra, dec));
CLUSTER raw_q3c_ang2pix_idx ON raw;
ANALYZE raw;