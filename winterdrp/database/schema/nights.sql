CREATE TABLE IF NOT EXISTS nights ( nid SERIAL, nightdate VARCHAR(20), nightid INT PRIMARY KEY,
    filename VARCHAR(500), filesize FLOAT,  checksum VARCHAR(50), status FLOAT);