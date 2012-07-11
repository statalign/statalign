package com.ppfold.algo;

import java.io.Serializable;

/**
 * Extended exponent datatype definition
 * 
 * @author Z.Sukosd
 */

public class PointRes implements Serializable {

	private static final long serialVersionUID = 6121110456253627605L;
	private static final double LOG_TWO = Math.log(2);
	// container for new number representation
	private float fraction; // = 0 OR >=1 0R <basis (in absolute value)
	private int exponent; // can be any integer, positive and negative

	public PointRes(float fraction, int exponent) {
		// Input fraction should always have exponent zero!!
		this.fraction = fraction;
		this.exponent = exponent;
	}
	
	public PointRes(PointRes point) {
		// Input fraction should always have exponent zero!!
		this.fraction = point.fraction;
		this.exponent = point.exponent;
	}

	public void setToFloat(float number) {
		this.exponent = 0;
		this.fraction = number;
		this.convert();
	}
	
	public void setToDouble(double number){
		long bits = Double.doubleToRawLongBits(number);
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		this.fraction = (float) fract;
		if (this.fraction == 0f) {
			this.fraction = 0;
			this.exponent = 0;
		} else {
			this.exponent = (int) ((bits >> 52) & 0x7FF) - 1023;
		}
	}

	public PointRes(double val) {
		long bits = Double.doubleToRawLongBits(val);
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		this.fraction = (float) fract;
		if (this.fraction == 0f) {
			this.fraction = 0;
			this.exponent = 0;
		} else {
			this.exponent = (int) ((bits >> 52) & 0x7FF) - 1023;
		}
	}

	@Override
	public PointRes clone() {
		PointRes newpoint = new PointRes(this.fraction, this.exponent);
		return newpoint.convert();
	}

	public void copyFrom(PointRes point) {
		this.exponent = point.exponent;
		this.fraction = point.fraction;
	}

	public int getExponent() {
		return exponent;
	}

	public void setExponent(int exponent) {
		this.exponent = exponent;
	}

	public float getFraction() {
		return fraction;
	}

	public void setFraction(float val) {
		this.fraction = val;
	}

	public PointRes convert() {
		if (fraction == 0f) {
			return this;
		}
		int bits = Float.floatToRawIntBits(fraction);
		int mantissaexp = ((bits >> 23) & 0xFF) - 127;
		exponent += mantissaexp;
		// set exp = 1 and then convert to double:
		fraction = Float.intBitsToFloat((bits & 0xBFFFFFFF) | 0x3f800000);
		return this;
	}

	public PointRes takeLogE(){
		double fraclog = Math.log(this.fraction);
		double exponentlog = LOG_TWO*this.exponent;
		this.fraction = (float) (fraclog + exponentlog);
		this.exponent = 0;
		this.convert();
		return this; 
	}
	
	public PointRes takeLog2(){
		double fraclog = Math.log(this.fraction)/LOG_TWO;
		this.fraction = (float) (fraclog + this.exponent);
		this.exponent = 0;
		this.convert();
		return this; 
	}
	
	public PointRes add(double d) {
		this.fraction = this.fraction + (float) d;
		this.convert();
		return this;
	}

	public PointRes add(PointRes point) {
		int diff = this.exponent - point.exponent;
		if ((diff < -126 && point.fraction != 0) || this.fraction == 0) {
			this.exponent = point.exponent;
			this.fraction = point.fraction;
		} else if ((diff > 127 && this.fraction != 0) || point.fraction == 0) {
			// do nothing, adding zero to the number
		} else if (point.exponent != this.exponent) {
			float newmantissa = point.fraction
					+ setExponent(this.fraction, diff);
			this.exponent = point.exponent;
			this.fraction = newmantissa;
			this.convert();
		} else {// equal exponents
			this.fraction = this.fraction + point.fraction;
			this.convert();
		}
		return this;
	}

	public boolean isSignificantlyLessThanZero(){
		if(this.fraction < 0f && this.exponent > -8){
			return true;
		}
		else return false; 
	}
	
	public boolean isLessThan(PointRes point) {
		if (this.fraction > 0f && point.fraction > 0f || this.fraction < 0f
				&& point.fraction < 0f) {
			if (this.exponent > point.exponent) {
				return false;
			} else if (this.exponent < point.exponent) {
				return true;
			} else { // case of equal exponents
				if (this.fraction < point.fraction) {
					return true;
				} else {
					return false;
				}
			}
		} else if (this.fraction >= 0f && point.fraction <= 0f) {
			return false;
		} else {
			return true;
		}
	}

	public boolean equals(PointRes point) {
		if (this.exponent == point.exponent && this.fraction == point.fraction) {
			return true;
		} else if (this.fraction == 0 && point.fraction == 0) {
			return true;
		} else {
			return false;
		}
	}

	public PointRes multiply(PointRes point) {
		this.fraction = this.fraction * point.fraction;
		this.exponent = this.exponent + point.exponent;
		this.convert();
		return this;
	}

	public PointRes multiply(double prob) {
		long bits = Double.doubleToRawLongBits(prob);
		// extract fraction part of probability
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		// cast to float
		float tmpfraction = (float) fract;
		if (tmpfraction == 0f) {
			// multiplying by zero
			this.fraction = 0;
			this.exponent = 0;
		} else {
			// multiplying by a non-zero value
			int tmpexponent = (int) ((bits >> 52) & 0x7FF) - 1023;
			this.fraction = this.fraction * tmpfraction;
			this.exponent = this.exponent + tmpexponent;
			this.convert();
		}
		return this;
	}

	public PointRes multiply(PointRes point, double prob) {
		this.multiply(point);
		this.multiply(prob);
		return this;
	}

	public PointRes divide(PointRes point) {
		if (point.fraction == 0f) {
			System.err.println("Division by zero!");
			return null;
		}
		this.fraction = this.fraction / point.fraction;
		this.exponent = this.exponent - point.exponent;
		this.convert();
		return this;
	}

	public float toFloat() {
		if (fraction == 0f) {
			return 0;
		} else if (-126 <= exponent && exponent <= 127) {
			// normal value
			return setExponent(fraction, exponent);
		} else if (exponent < -126) {
			return 0;
		} else {
			return fraction > 0f ? Float.POSITIVE_INFINITY
					: Float.NEGATIVE_INFINITY;
		}

	}

	public double toDouble() {
		if (fraction == 0f) {
			return 0;
		} else if (-126 <= exponent && exponent <= 127) {
			// normal value
			return (double) setExponent(fraction, exponent);
		} else if (exponent < -126) {
			return 0;
		} else {
			return fraction > 0f ? Double.POSITIVE_INFINITY
					: Double.NEGATIVE_INFINITY;
		}
	}

	public String toString() {
		return "" + this.fraction + " x 2^" + this.exponent;
	}

	// exponent must be between -126 and 127
	private static float setExponent(float d, int expin) {
		int bits = Float.floatToRawIntBits(d);
		bits = bits & 0x807fffff; // exp bits are 0
		bits = bits | ((expin + 127) << 23);
		return Float.intBitsToFloat(bits);
	}

	public void subtract(PointRes pr, PointRes tmp) {
		tmp.copyFrom(pr);
		tmp.multiply(-1);
		this.add(tmp);
	}
}