package statalign.utils;

/**
 * Provides simple statistical information such as average and standard deviation for a set of values.
 * These are calculated dynamically without storing the data values themselves. The applied formulas
 * are optimized for numerical stability.
 *
 * @author novak.adam
 *
 */
public class SimpleStats {

	private String name;

	private int n;
	private double avg;			// = 1/n*sum^n{x_i}
	private double dsqSum;		// = sum^n{(x_i-avg)^2}

	private double min = Double.MAX_VALUE;
	private double max = Double.MIN_VALUE;

	public SimpleStats() {
	}

	public SimpleStats(String name) {
		this.name = name;
	}

	public void addData(double value) {
		double d = value-avg;
		n++;
		avg += d/n;
		dsqSum += (n-1)*d*d/n;
		if(value < min) {
			min = value;
		}
		if(value > max) {
			max = value;
		}
	}

	public void addData(double[] values) {
		for(double value : values) {
			addData(value);
		}
	}

	public void addData(int[] values) {
		for(int value : values) {
			addData(value);
		}
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getN() {
		return n;
	}

	public double getAvg() {
		return avg;
	}

	public double getStdDev() {
		if(n <= 1) {	// to avoid division by zero
			return 0;
		}
		return Math.sqrt(dsqSum/(n-1));
	}

	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	public void reset() {
		n = 0;
		avg = 0;
		dsqSum = 0;
	}

	@Override
	public String toString() {
		String res = "";

		if(name != null) {
			res = name+": ";
		}

		return res+"n="+getN()+" avg="+getAvg()+" stddev="+getStdDev()+" min="+min+" max="+max;
	}

	public static void main(String[] args) {
		SimpleStats test = new SimpleStats("test");
		test.addData(.5);
		test.addData(.10);
		test.addData(.7);
		System.out.println(test);
		test.reset();
		System.out.println(test);
		test.addData(new double[] {4,5});
		System.out.println(test);
	}
}
