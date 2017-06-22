function cont = KEGGReactioncontainsGlycan(equation)
    cont = false;
    match = regexp(equation, 'G\d{5}', 'match');
    if isempty(match)
        cont = false;
        return
    else
        cont = true;
    end
end