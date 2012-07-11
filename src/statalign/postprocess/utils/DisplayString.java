package statalign.postprocess.utils;
/**
 * A class that stores the strings to be displayed and thus enable zooming on the consensus trees / networks
 * 
 * @author wood
 *
 */
public class DisplayString {
	public String label;
	public int x;
	public int y;
	public DisplayString(String label,int x, int y){
		this.x = x;
		this.y = y;
		this.label = label;
	}

}
