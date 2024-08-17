% estimating distance spectrum for turbo code

distance_data = importdata("../output/rank1000_metrics.txt");

h = histogram(distance_data/4);

binCounts  = h.Values;
binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end)) / 2;

for i = 1:length(binCounts)
    text(binCenters(i), binCounts(i), num2str(binCounts(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

xticks([0, 5, 10, 13, 14, 15, 20, 24]);
xlim([0, 24]);
xlabel('Hamming distance', 'FontSize', 15);
ylabel("Count", 'FontSize', 15);
title('Estimation of distance spectrum', 'FontSize', 17);


%%

metric_diff = [30.5022, 12.7782, 12.4392, 28.1851, 23.1182, 22.3352, 13.8581, 32.4074, 30.3837, 17.18, 28.532, 12.1856, 23.7634, 20.8609, 23.3883, 24.0023, 28.2641, 22.6156, 11.5301, 26.5709, 26.4288, 25.3233, 30.6608, 33.9848, 20.0199, 17.5467, 14.7557, 8.25915, 7.91673, 35.2296, 29.7888, 27.9521, 32.7118, 29.0944, 30.2164, 22.7089, 23.102, 32.8893, 24.5674, 23.7138, 27.6124, 17.8881, 27.318, 17.2988, 37.2012, 11.5131, 24.5545, 16.4131, 24.6216, 29.717, 18.9551, 22.2041, 42.3978, 34.8277, 20.9201, 30.278, 34.2611, 34.1605, 10.2256, 27.3826, 30.7636, 9.7204, 16.899, 38.2259, 22.4344, 20.2869, 21.9691, 33.0745, 17.9103, 18.597, 24.2984, 16.7152, 34.1272, 19.9136, 23.7253, 23.6912, 32.5275, 39.3548, 12.6315, 42.1628, 20.0149, 40.4615, 20.1846, 33.894, 40.9309, 25.4197, 29.966, 27.2394, 29.5217, 18.1901, 20.9558, 10.9823, 48.5857, 18.2109, 38.0195, 25.8564, 31.7778, 27.4771, 45.9177, 17.2782];

histogram(metric_diff);
xlabel("Return metric - true metric", 'FontSize', 15);
ylabel("Count", 'FontSize', 15);
title('Error Event Metric difference', 'FontSize', 17);
