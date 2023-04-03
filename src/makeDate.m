function [y] = makeDate(x)
    y = datetime(x, 'ConvertFrom', 'posixtime', 'Format', ...
        'yyyy/MM/dd (DDD) HH:mm', 'TimeZone', 'America/Los_Angeles');
end

