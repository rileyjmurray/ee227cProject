% world coloring tasks for each interesting set of continents
% color of x(1) needs to be different from color x(2).
relation = @(x) x(1) ~= x(2);

numColors = 4; %change this

% AMERICAS
% 1. alaska
% 2. canada
% 3. usa
% 4. mexico
% 5. guatemala
% 6. belize
% 7. el salvador
% 8. honduras
% 9. nicaragua
% 10. costa rica
% 11. panama
% 12. colombia
% 13. venezuela
% 14. guyana
% 15. suriname
% 16. French guyana 
% 17. Brazil
% 18. Ecuador
% 19. Peru
% 20. Bolivia
% 21. Paraguay
% 22. Chile
% 23. Argentina
% 24. Uruguay

borders = ...
    [1 2;...
    2 3;...
    3 4;...
    4 5;...
    4 6;...
    5 6;...
    5 7;...
    5 8;...
    7 8;...
    8 9;...
    9 10;...
    10 11;...
    11 12;...
    12 13;...
    12 17;...
    12 18;...
    12 19;...
    13 14;...
    13 17;...
    14 15;...
    14 17;...
    15 16;...
    15 17;...
    16 17;...
    17 19;...
    17 20;...
    17 21;...
    17 23;...
    17 24;...
    18 19;...
    19 20;...
    19 22;...
    20 21;...
    20 22;...
    20 23;...
    21 23;...
    22 23;...
    23 24];

constraints = cell( size( borders, 1 ), 1 );
for i = 1 : size( borders, 1 )
    scope = borders( i, : );
    constraints{ i } = Constraint( scope, relation );
end

americasProblem = CSP(ones(length(borders),1), constraints, 0:(numColors-1), 24);
clear constraints i scope;
%%
% AFRICA
% 1. Morocco
% 2. Algeria
% 3. Tunesia
% 4. Libia
% 5. Egypt
% 6. Mauritania
% 7. Mali
% 8. Niger
% 9. Chad
% 10. Sudan
% 11. Eritrea
% 12. Senegal
% 13. Gambia
% 14. Guinea Bissau
% 15. Guinea
% 16. Sierra Leona
% 17. Liberia
% 18. Cote d'Ivoire
% 19. Burkina Faso
% 20. Ghana
% 21. Togo
% 22. Benin
% 23. Nigeria
% 24. Cameroon
% 25. South Sudan
% 26. Ethiopia
% 27. Djibouti
% 28. Somalia
% 29. Equatorial Guinea
% 30. Gabon
% 31. CAR
% 32. Congo
% 33. DR Congo
% 34. Uganda
% 35. Rwanda
% 36. Burundi
% 37. Kenya
% 38. Tanzania
% 39. Angola
% 40. Zambia
% 41. Malawi
% 42. Mozambique
% 43. Madagascar
% 44. Namibia
% 45. Botswana
% 46. Zimbabwe
% 47. South Africa
% 48. SwaziLand
% 49. Lesotho

borders = [
    1 2; 1 6; 2 3; 2 4; 2 6; 2 7; 2 8; 3 4; 4 5; 4 8; 4 9; 4 10;
    5 10; 6 12; 6 7; 7 8; 7 12; 7 15; 7 18; 7 19; 8 9; 8 19; 8 22; 8 23;
    12 13; 12 14; 12 15; 14 15; 15 16; 15 17; 15 18; 16 17; 17 18; 
    18 19; 18 20; 19 20; 19 21; 19 22; 20 21; 21 22; 22 23;
    9 10; 9 23; 9 24; 9 25; 9 31; 10 31; 10 11; 10 25; 10 26;
    11 26; 11 27; 23 24; 24 31; 24 29; 24 30; 24 32; 25 26; 25 33;
    25 34; 25 37; 26 37; 26 27; 26 28; 27 28; 28 37;
    29 30; 30 32; 32 33; 33 34; 33 35; 33 36; 33 38; 33 39; 33 40;
    34 37; 34 35; 34 38; 35 36; 35 38; 36 38; 37 38; 38 40; 38 41;
    38 42; 39 40; 39 44; 40 42; 40 46; 40 41; 40 44; 41 42;
    42 43; 42 48; 42 47; 42 26; 44 45; 44 47; 44 46; 45 46; 45 47;
    46 47; 47 48; 47 49 ];

constraints = cell( size( borders, 1 ), 1 );

for i = 1 : size( borders, 1 )
    scope = borders( i, : );
    constraints{ i } = Constraint( scope, relation );
end

africaProblem = CSP( ones( size( borders, 1 ), 1 ), constraints, 0:(numColors-1), 49);
clear relation;
clear constraints i scope;