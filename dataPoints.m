function points = dataPoints(Vc, margin)
    part1 = linspace(0.01, Vc-margin, 25);
    part2 = linspace(Vc-margin, Vc+margin, 100);
    part3 = linspace(Vc+margin, 2*Vc, 25);
    points = [part1, part2, part3];
end

