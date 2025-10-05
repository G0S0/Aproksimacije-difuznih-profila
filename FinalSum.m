ApproximationBlue
ApproximationGreen
ApproximationRed

sumOfRGBGaussians = @(r)normalizedGaussBlue(r) + normalizedGaussRed(r) + normalizedGaussGreen(r);
fplot(@(rs) sumOfRGBGaussians(rs), [0, 4], 'Color', [0, 1, 1], 'LineWidth', 2) %cyan


sumOfRGBGaussians(1.0)
sumOfRGBGaussians(0.75)
sumOfRGBGaussians(0.5)
sumOfRGBGaussians(0.25)
sumOfRGBGaussians(0.5)
sumOfRGBGaussians(0.75)
sumOfRGBGaussians(1.0)

sumOfRGBGaussianWeights = (sumOfRGBGaussians(0.25) + 2* sumOfRGBGaussians(0.5) + ...
    2* sumOfRGBGaussians(0.75) + 2*sumOfRGBGaussians(1.0))

weight025 = sumOfRGBGaussians(0.25)/sumOfRGBGaussianWeights
weight050 = sumOfRGBGaussians(0.50)/sumOfRGBGaussianWeights
weight075 = sumOfRGBGaussians(0.75)/sumOfRGBGaussianWeights
weight100 = sumOfRGBGaussians(1)/sumOfRGBGaussianWeights

sumOfWeights = 2*weight100 + 2*weight075 + weight025 + 2*weight050