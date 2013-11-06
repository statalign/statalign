package statalign.postprocess;

import java.awt.Color;

public class Track {
	public Color color;
	public double[] scores;
	public Track(Color _color, double[] _scores) {
		color = _color;
		scores = _scores;
	}
}
