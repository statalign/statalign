package statalign.postprocess;

import java.awt.Color;

public class Track {
	public Color color;
	public double[] scores;
	public double min;
	public double max;
	public double mean;
	public Track(Color _color, double[] _scores) {
		color = _color;
		scores = _scores;
	}
}
