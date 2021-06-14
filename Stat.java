import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class Stat {

	public static double Mean(double[] input) {
		int sum = 0;
		for (double a : input) {
			sum += a;
		}
		return sum / input.length;
	}

	public static double StandardDeviation(double[] input) {
		double Mean = Mean(input);
		int sum = 0;
		for (double a : input) {
			sum += Math.pow((a - Mean), 2);
		}
		return Math.sqrt(sum / input.length);
	}

	public static double Variation(double[] input) {
		double Mean = Mean(input);
		int sum = 0;
		for (double a : input) {
			sum += Math.pow((a - Mean), 2);
		}
		return sum / (input.length - 1);
	}

	public static double SEOfMean(double[] input) {
		return StandardDeviation(input) / (Math.sqrt(input.length));

	}

	public static double CoefficientOfVariation(double[] input) {
		return (StandardDeviation(input) / Mean(input)) * 100;
	}

	public static double FirstQuartile(double[] input) {
		double[] temp = input.clone();
		Arrays.sort(temp);
		int length = temp.length + 1;
		double position = length * (.25);
		int lowerBound = ((int) position) - 1;
		int upperBound = ((int) position);
		if (lowerBound == length * (.25)) {
			return temp[lowerBound];
		}
		double difference = (temp[upperBound] - temp[lowerBound]) * (position - (lowerBound + 1));
		return temp[lowerBound] + difference;
	}

	public static double Median(double[] input) {
		double[] temp = input.clone();
		Arrays.sort(temp);
		int length = temp.length + 1;
		double position = length * (.5);
		int lowerBound = ((int) position) - 1;
		int upperBound = ((int) position);
		if (lowerBound == length * (.5)) {
			return temp[lowerBound];
		}
		double difference = (temp[upperBound] - temp[lowerBound]) * (position - (lowerBound + 1));
		return temp[lowerBound] + difference;
	}

	public static double ThirdQuartile(double[] input) {
		double[] temp = input.clone();
		Arrays.sort(temp);
		int length = temp.length + 1;
		double position = length * (.75);
		int lowerBound = ((int) position) - 1;
		int upperBound = ((int) position);
		if (lowerBound == length * (.75)) {
			return temp[lowerBound];
		}
		double difference = (temp[upperBound] - temp[lowerBound]) * (position - (lowerBound + 1));
		return temp[lowerBound] + difference;
	}

	public static double InterQuartileRange(double[] input) {
		return ThirdQuartile(input) - FirstQuartile(input);
	}

	public static double[] Mode(double[] input) {
		HashMap<Double, Double> Frequency = new HashMap<Double, Double>();
		for (double a : input) {
			if (!Frequency.containsKey(a)) {
				Frequency.put(a, 1.0);
			} else {
				Frequency.replace(a, Frequency.get(a) + 1.0);
			}
		}
		double largestNumber = 0;
		int numOfLarge = 0;
		for (Map.Entry<Double, Double> entry : Frequency.entrySet()) {
			double value = entry.getValue();
			if (value > largestNumber) {
				largestNumber = value;
				numOfLarge = 1;
			} else if (value == largestNumber) {
				numOfLarge++;
			}
		}
		double[] output = new double[numOfLarge];
		int position = 0;
		for (Map.Entry<Double, Double> entry : Frequency.entrySet()) {
			double key = entry.getKey();
			double value = entry.getValue();
			if (value == largestNumber) {
				output[position] = key;
				position++;
			}
		}
		Arrays.sort(output);
		return output;
	}

	public static double Minimum(double[] input) {
		double min = input[0];
		for (double a : input) {
			if (a < min) {
				min = a;
			}
		}
		return min;
	}

	public static double Maximum(double[] input) {
		double max = input[0];
		for (double a : input) {
			if (a > max) {
				max = a;
			}
		}
		return max;
	}

	public static double Range(double[] input) {
		return Maximum(input) - Minimum(input);
	}

	public static double Factorial(int input) {
		double output = 1.0;
		if (input == 0) {
			return 1.0;
		}
		while (input >= 1) {
			output = input * output;
			input--;
		}
		return output;
	}

	// Cecil Hastings, Jr., Approximations for Digital Computers, Princeton, NJ:
	// Princeton University Press 1955, pp. 187
	public static double ErrorFunction(double input) {
		return 1 - 1 / Math.pow((1 + .0705230784 * input + .0422820123 * Math.pow(input, 2)
				+ .0092705272 * Math.pow(input, 3) + .0001520143 * Math.pow(input, 4) + .0002765672 * Math.pow(input, 5)
				+ .0000430638 * Math.pow(input, 6)), 16);
	}

	public static double ZScoreToPercentile(double input, boolean twoTailed, boolean rightTailedOrOuterTailed) {
		double error = ErrorFunction(input / Math.sqrt(2));
		double percentile = error / 2.0 + .5;

		if (twoTailed == false) {
			if (rightTailedOrOuterTailed == true) {
				return 1.0 - percentile;
			} else {
				return percentile;
			}
		} else {
			if (rightTailedOrOuterTailed == true) {
				if (percentile > .5) {
					return (1.0 - percentile) * 2.0;
				} else {
					return percentile * 2.0;
				}
			} else {
				if (percentile > .5) {
					return (percentile - .5) * 2.0;
				} else {
					return (.5 - percentile) * 2.0;
				}
			}
		}
	}

	/*
	 * Approximation of Inverse Error Function by Peter John Acklam Modified code of
	 * javascript implementation by Alankar Misra
	 * https://web.archive.org/web/20150915095009/http://home.online.no/~pjacklam/
	 * notes/invnorm/impl/misra/normsinv.html
	 */
	public static double ProbitFunction(double p) {
		// Coefficients in rational approximations
		double[] a = new double[] { -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
				1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00 };

		double[] b = new double[] { -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
				6.680131188771972e+01, -1.328068155288572e+01 };

		double[] c = new double[] { -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
				-2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };

		double[] d = new double[] { 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
				3.754408661907416e+00 };

		// Define break-points.
		double plow = 0.02425;
		double phigh = 1 - plow;

		// Rational approximation for lower region:
		if (p < plow) {
			double q = Math.sqrt(-2 * Math.log(p));
			return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
					/ ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
		}

		// Rational approximation for upper region:
		if (phigh < p) {
			double q = Math.sqrt(-2 * Math.log(1 - p));
			return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
					/ ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
		}

		// Rational approximation for central region:
		double q = p - 0.5;
		double r = q * q;
		return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
				/ (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
	}

	public static double PercentileToZscore(double percentile, boolean twoTailed, boolean rightTailedOrOuterTailed) {

		if (twoTailed == false) {
			if (rightTailedOrOuterTailed == true) {
				return ProbitFunction(1.0 - percentile);
			} else {
				return ProbitFunction(percentile);
			}
		} else {
			if (rightTailedOrOuterTailed == true) {
				return ProbitFunction(1.0 - (percentile / 2.0));
			} else {
				return ProbitFunction(percentile / 2.0 + .5);
			}
		}

	}

	static double GammaFunction(double x) {
		double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1) + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
				+ 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
		return Math.exp(tmp + Math.log(ser * Math.sqrt(2 * Math.PI)));
	}

	static double TDistributionPDF(double tScore, double degreesOfFreedom) {
		return (GammaFunction((degreesOfFreedom + 1.0) / 2.0)
				/ (Math.sqrt(degreesOfFreedom * Math.PI) * GammaFunction(degreesOfFreedom / 2.0)))
				* Math.pow(1.0 + (Math.pow(tScore, 2.0) / degreesOfFreedom), -(degreesOfFreedom + 1.0) / 2.0);
	}

	static double TScoreToPerecentile(double tScore, double degreesOfFreedom, boolean twoTailed,
			boolean rightTailedOrOuterTailed) {
		boolean left = false;
		if (tScore < 0) {
			tScore *= -1.0;
			left = true;
		}

		DecimalFormat df = new DecimalFormat("#.#####");
		tScore = Double.parseDouble(df.format(tScore));
		double percentile = 0.0;
		for (double a = tScore; a > 0; a -= .00001) {
			percentile += .00001 * TDistributionPDF(a, degreesOfFreedom);
		}
		percentile += .5;
		if (percentile > 1.0) {
			percentile = 1.0;
		}
		if (left) {
			percentile = 1.0 - percentile;
		}
		if (tScore == 0.0) {
			percentile = .5;
		}

		if (twoTailed == false) {
			if (rightTailedOrOuterTailed == true) {
				return 1.0 - percentile;
			} else {
				return percentile;
			}
		} else {
			if (rightTailedOrOuterTailed == true) {
				if (percentile > .5) {
					return (1.0 - percentile) * 2.0;
				} else {
					return percentile * 2.0;
				}
			} else {
				if (percentile > .5) {
					return (percentile - .5) * 2.0;
				} else {
					return (.5 - percentile) * 2.0;
				}
			}
		}
	}

	static double PercentileToTScoreHelper(double percentile, double degreesOfFreedom) {
		double perc = percentile;	
		double ZScore = PercentileToZscore(perc,false,false);
		boolean forwards = true;
		double xIncrement = 1.0;
		if(percentile < .5) {
			forwards = false;
			xIncrement = -1.0;
		}
		double xPosition = ZScore;
		boolean continueLoops = true;
		double currentPerc = TScoreToPerecentile(xPosition, degreesOfFreedom, false, false);
		while(continueLoops) {
			if(perc-.00001 <= currentPerc && perc+.00001 >= currentPerc ) {
				continueLoops = false;
			}
			else {
				if(forwards) {
					if(currentPerc > perc) {
						forwards = false;
						xIncrement *= -0.1;
					}
					else {
						xPosition += xIncrement;
						currentPerc = TScoreToPerecentile(xPosition, degreesOfFreedom, false, false);
					}
				}
				else {
					if(currentPerc < perc) {
						forwards = true;
						xIncrement *= -0.1;
					}
					else {
						xPosition += xIncrement;
						currentPerc = TScoreToPerecentile(xPosition, degreesOfFreedom, false, false);
					}
				}
			}
		}
		return xPosition;
	}

	static double PercentileToTScore(double percentile, double degreesOfFreedom, boolean twoTailed, boolean rightTailedOrOuterTailed) {
		if(twoTailed == false) {
			if(rightTailedOrOuterTailed == true) {
				return PercentileToTScoreHelper(1.0-percentile, degreesOfFreedom);
			}
			else {
				return PercentileToTScoreHelper(percentile, degreesOfFreedom);
			}
		}
		else {
			if(rightTailedOrOuterTailed == true) {
				return PercentileToTScoreHelper(1.0 - (percentile / 2.0), degreesOfFreedom);
			}
			else {
				return PercentileToTScoreHelper(percentile / 2.0 + .5, degreesOfFreedom);
			}
		}	
	}
	public static double OneSampleZHypothesis(double mean, double stdev, double hypothesis, boolean twoTailed,
			Boolean rightTailed) {
		return ZScoreToPercentile((hypothesis - mean) / stdev, twoTailed, rightTailed);
	}

	public static double TwoSamplePearsonsCorrelationConstant(double[] setX, double[] setY) {
		double top = 0;
		double SummationOne = 0;
		for (int a = 0; a < setX.length; a++) {
			SummationOne += setX[a] * setY[a];
		}
		double SummationTwo = 0;
		for (int a = 0; a < setX.length; a++) {
			SummationTwo += setX[a];
		}
		double SummationThree = 0;
		for (int a = 0; a < setX.length; a++) {
			SummationThree += setY[a];
		}
		top = setX.length * SummationOne - SummationTwo * SummationThree;
		double bottom = 0;
		double SummationFour = 0;
		for (int a = 0; a < setX.length; a++) {
			SummationFour += Math.pow(setX[a], 2);
		}
		double SummationFive = 0;
		for (int a = 0; a < setX.length; a++) {
			SummationFive += Math.pow(setY[a], 2);
		}
		bottom = Math.sqrt((setX.length * SummationFour - Math.pow(SummationTwo, 2))
				* (setX.length * SummationFive - Math.pow(SummationThree, 2)));
		return top / bottom;
	}

	public static double TwoSampleRegressionSlope(double[] setX, double[] setY, double r) {
		return r * (StandardDeviation(setY) / StandardDeviation(setX));
	}

	public static double TwoSampleRegressionIntercep(double[] setX, double[] setY, double slope) {
		return Mean(setY) - slope * Mean(setX);
	}

	public static double TwoSampleResidualStandardDeviation(double[] setX, double[] setY, double intercep,
			double slope) {
		double top = 0;
		for (int a = 0; a < setX.length; a++) {
			top += Math.pow(setY[a] - intercep + slope * setX[a], 2);
		}
		return Math.sqrt(top / (setX.length - 2));
	}

	public static double TwoSampleStandardErrorOfTheSlope(double[] setX, double ResidualStandardDeviation) {
		return ResidualStandardDeviation / (Math.sqrt(setX.length - 1) * StandardDeviation(setX));
	}

	public static double TwoSampleStandardErrorOfTheIntercep(double[] setX, double ResidualStandardDeviation) {
		double summation = 0;
		double meanX = Mean(setX);
		for (int a = 0; a < setX.length; a++) {
			summation += Math.pow(setX[a] - meanX, 2);
		}
		return ResidualStandardDeviation * Math.sqrt(1 / setX.length + Math.pow(meanX, 2) / summation);
	}

	public static double TwoSampleRegressionPredictedMeanValueStandardError(double x, double[] setX, double ResidualSD,
			double SESlope) {
		return Math.sqrt(Math.pow(SESlope, 2) * Math.pow(x - Mean(setX), 2) + (Math.pow(ResidualSD, 2) / setX.length));
	}

	public static double TwoSampleRegressionPredictedIndividualValueStandardError(double x, double[] setX,
			double ResidualSD, double SESlope) {
		return Math.sqrt(Math.pow(SESlope, 2) * Math.pow(x - Mean(setX), 2) + (Math.pow(ResidualSD, 2) / setX.length)
				+ Math.pow(ResidualSD, 2));
	}

	public static double TwoSampleRegressionPredictY(double x, double intercep, double slope) {
		return intercep + slope * x;
	}
	
	public static double SEOfTheMean(double StDev, double size) {
		return StDev/Math.sqrt(size);
	}
	
	public static double SEOfDifferenceOfTwoMeans(double stDevOne, double sizeOne, double stDevTwo, double sizeTwo) {
		return Math.sqrt(Math.pow(stDevOne,2)/sizeOne+Math.pow(stDevTwo,2)/sizeTwo);
	}
	
	public static double DFOfDifferenceOfTwoMeans(double stDevOne, double sizeOne, double stDevTwo, double sizeTwo) {
		return Math.pow(Math.pow(stDevOne, 2)/sizeOne + Math.pow(stDevTwo, 2)/sizeTwo,2)/
				((1.0/(sizeOne-1.0))*(Math.pow(stDevOne, 2)/sizeOne)+(1.0/(sizeTwo-1.0))*(Math.pow(stDevTwo, 2)/sizeTwo));
	}
	
	public static double pooledStandardDeviationOfTwoMeans(double stDevOne, double sizeOne, double stDevTwo, double sizeTwo) {
		return ((sizeOne-1.0)*Math.pow(stDevOne,2) + (sizeTwo-1.0)*Math.pow(stDevTwo,2))/(sizeOne+sizeTwo-2.0);
	}
	
	public static double pooledSEOfTwoMeans(double pooledStDev, double sizeOne, double sizeTwo) {
		return pooledStDev*Math.sqrt(1.0/sizeOne + 1.0/sizeTwo);
	}
	
	public static double DFOfPooledMeans(double sizeOne, double sizeTwo) {
		return sizeOne+sizeTwo-2.0;
	}
	
	
	public static double SEOfProportion(double proportion, double size) {
		return Math.sqrt((proportion*(1.0-proportion))/size);
	}
	
	public static double SEOfDifferenceOfTwoProportions(double proportionOne, double sizeOne, double proportionTwo, double sizeTwo) {
		return Math.sqrt((proportionOne*(1.0-proportionOne))/sizeOne + (proportionTwo*(1.0-proportionTwo))/sizeTwo);
	}
	
	public static double proportionPooled(double proportionOne, double sizeOne, double proportionTwo, double sizeTwo) {
		return (proportionOne+proportionTwo)/(sizeOne+sizeTwo);
	}
	
	public static double SEOfPooledProportions(double pooledProportion, double sizeOne, double sizeTwo) {
		return Math.sqrt((pooledProportion*(1.0-pooledProportion))/sizeOne + (pooledProportion*(1.0-pooledProportion))/sizeTwo);
	}
}
