function values = parse_integer_list(text, lowerBound, upperBound, allowEmpty)
%PARSE_INTEGER_LIST Parse comma/space lists, with optional outer brackets.
if nargin < 4, allowEmpty = false; end
text = strtrim(char(text));
if startsWith(text,'[') && endsWith(text,']')
    text = strtrim(text(2:end-1));
end
if isempty(text) && allowEmpty
    values = [];
    return
end
if isempty(text) || isempty(regexp(text,'^\d+(?:[\s,]+\d+)*$','once'))
    error('edgepreserve:integerList', 'Use integer IDs separated by commas or spaces.');
end
values = str2double(regexp(text,'\d+','match'));
if any(~isfinite(values) | values < lowerBound | values > upperBound)
    error('edgepreserve:integerRange', 'Values must be integers in [%g,%g].', ...
        lowerBound,upperBound);
end
values = unique(values);
end
