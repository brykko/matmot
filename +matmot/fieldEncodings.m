function encodings = fieldEncodings()
%FIELDENCODINGS

FIELDS = {
    'frameIdx'
    'frameTimestamp'
    'frameLatency'
    'pos'
    'rot'
    'posError'
    'posTracked'};

BYTE_INDS = {
    1:4
    5:12
    13:16
    17:28
    29:44
    45:48
    49};

ENCODINGS = {
    'int32'
    'double'
    'single'
    'single'
    'single'
    'single'
    'uint8'};

encodings = struct('field', FIELDS, 'byte_inds', BYTE_INDS, 'encoding', ENCODINGS);

end

