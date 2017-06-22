function [id, rp] = parseRPAIR(rpair)

tmp = textscan(rpair, '%s', 'delimiter', '\n', 'whitespace', '');

% find RDM pattern
entry = strfind(tmp{1}, 'RDM');
ind = find(~cellfun('isempty', entry), 1);
rdm = tmp{1}(ind+1);
token ='1';
a = regexprep(rdm,token,'');
rdm = strtrim(a);
rdmpattern = textscan(rdm{1}, '%s%s%s', 'delimiter', ':', 'whitespace', '');

% find type
entry = strfind(tmp{1}, 'TYPE');
ind = find(~cellfun('isempty', entry), 1);
type = strtrim(tmp{1}(ind)); 
token ='^TYPE';
a = regexprep(type,token,'');
type = strtrim(a);
type = type{:};
type = regexp(type,'\s','split');

% find compounds
nameid = strfind(tmp{1}, 'NAME');
ind = ~cellfun('isempty', nameid);
n = strtrim(tmp{1}(ind));
name = strtrim(regexprep(n,'NAME',''));
compounds = regexp(name,'C\d*', 'match');
compounds = compounds{:};%todo: names needed?

%find reactions
reactions = regexp(rpair,'R\d{5}', 'match');

%find id
index = strfind(tmp{1}, 'ENTRY');
ind = find(~cellfun('isempty', index), 1);
entry = strtrim(tmp{1}(ind));
id = regexp(entry,'RP\d{5}', 'match');
id = id{:};
id = id{:};

%compile data 
name = name{:};
rp.id = id;
rp.name = name;
rp.compounds = compounds;
rp.reactions = reactions;
rp.type = type;
rp.rdmpattern = rdmpattern;


end  