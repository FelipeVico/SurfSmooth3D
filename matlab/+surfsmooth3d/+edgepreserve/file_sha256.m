function hash = file_sha256(filename)
%FILE_SHA256 Hash input bytes in bounded blocks without modifying the file.
fid = fopen(filename,'rb');
if fid < 0, error('edgepreserve:hashOpen', 'Cannot read %s.',filename); end
cleanup = onCleanup(@() fclose(fid));
digest = java.security.MessageDigest.getInstance('SHA-256');
while true
    block = fread(fid,1024*1024,'*uint8');
    if isempty(block), break; end
    digest.update(typecast(block,'int8'));
end
bytes = typecast(digest.digest(),'uint8');
hash = lower(reshape(dec2hex(bytes,2).',1,[]));
end
