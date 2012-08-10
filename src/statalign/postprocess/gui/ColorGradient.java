package statalign.postprocess.gui;

import java.awt.Color;
import java.util.Arrays;

/**
 *
 * @author Michael Golden
 */
public class ColorGradient {

    public Color[] colours;
    public float[] positions;

    public ColorGradient(Color c1, Color c2) {
        colours = new Color[2];
        colours[0] = c1;
        colours[1] = c2;

        positions = new float[2];
        positions[0] = 0;
        positions[1] = 1;
    }

    public ColorGradient(Color c1, Color c2, Color c3) {
        colours = new Color[3];
        colours[0] = c1;
        colours[1] = c2;
        colours[2] = c3;

        positions = new float[3];
        positions[0] = 0f;
        positions[1] = 0.5f;
        positions[2] = 1f;
    }

    public ColorGradient(Color[] colours) {
        this.colours = colours;

        positions = new float[colours.length];
        for (int i = 0; i < colours.length; i++) {
            positions[i] = i / ((float) (colours.length - 1));
        }

    }

    public void distributeColors() {
        for (int i = 0; i < colours.length; i++) {
            positions[i] = i / ((float) (colours.length - 1));
        }
    }

    public void reverseOrder()
    {
        Color [] tempColors = new Color[colours.length];
        for(int i = 0 ; i < tempColors.length ; i++)
        {
            tempColors[colours.length-1-i] = colours[i];
        }
        float [] tempPositions = new float[positions.length];
        for(int i = 0 ; i < tempPositions.length ; i++)
        {
            tempPositions[positions.length-1-i] = 1-positions[i];
        }
        colours = tempColors;
        positions = tempPositions;
    }

    public ColorGradient(Color[] colours, float[] positions) {
        this.colours = colours;
        this.positions = positions;
    }

    public Color getColor(float val) {
        if (val > 1) {
            return colours[colours.length - 1];
        } else if (val < 0) {
            return colours[0];
        }

        int c1 = 0;
        int c2 = colours.length - 1;
        float lower = 0;
        float upper = 1;
        int pos = 0;
        for (int i = 1; i < positions.length; i++) {
            lower = positions[i - 1];
            upper = positions[i];
            if (lower <= val && val <= upper) {
                pos = i;
                break;
            }
        }

        c1 = Math.max(pos - 1, 0);
        c2 = pos;
        if (val > positions[positions.length - 1]) {
            c1 = positions.length - 1;
            c2 = positions.length - 1;
        }


        float perc = (val - lower) / (upper - lower);
        int r = colours[c1].getRed() + ((int) (perc * (colours[c2].getRed() - colours[c1].getRed())));
        int g = colours[c1].getGreen() + ((int) (perc * (colours[c2].getGreen() - colours[c1].getGreen())));
        int b = colours[c1].getBlue() + ((int) (perc * (colours[c2].getBlue() - colours[c1].getBlue())));
        int a = colours[c1].getAlpha() + +((int) (perc * (colours[c2].getAlpha() - colours[c1].getAlpha())));

        return new Color(r, g, b, a);
    }

    public String toString() {
        String ret = "";
        for (int i = 0; i < colours.length; i++) {
            ret += "rgba(" + colours[i].getRed() + "," + colours[i].getGreen() + "," + colours[i].getBlue() + "," + colours[i].getAlpha() + "):"+positions[i];
            ret += ";";
        }
        return ret;

    }

    public static ColorGradient getValue(String fromString) {
        String[] colorSplit = fromString.split(";");
        Color[] colors = new Color[colorSplit.length];
        float[] positions = new float[colorSplit.length];
        for (int i = 0; i < colors.length; i++) {
            String [] colorSplit2 = colorSplit[i].split(":");
            if(colorSplit2.length >= 1)
            {
                String[] rgba = colorSplit2[0].replaceAll("([^0-9,])+", "").split(",");
                colors[i] = new Color(Integer.parseInt(rgba[0]), Integer.parseInt(rgba[1]), Integer.parseInt(rgba[2]), Integer.parseInt(rgba[3]));
            }
            if(colorSplit2.length >= 2)
            {
                positions[i] = Float.parseFloat(colorSplit2[1]);
            }
            else
            {
                positions[i] = i / ((float) (colors.length - 1));
            }
        }
        return new ColorGradient(colors, positions);
    }

    public ColorGradient clone ()
    {
        return new ColorGradient(Arrays.copyOf(colours, colours.length), Arrays.copyOf(positions, positions.length));
    }

    public static void main(String[] args) {
        Color[] colours = {Color.green, Color.white, Color.blue};
        System.out.println((new ColorGradient(colours)).getColor(0.99f));
    }
}
