package fr.orsay.lri.varna.models.annotations;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.controlers.ControleurClicMovement;
import fr.orsay.lri.varna.models.VARNAConfigLoader;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.RNA;
import fr.orsay.lri.varna.views.VueHighlightRegionEdit;
import fr.orsay.lri.varna.views.VueUI;

public class HighlightRegionAnnotation implements Serializable {
	private static final long serialVersionUID = 7087014168028684775L;
	public static final Color DEFAULT_OUTLINE_COLOR = Color.decode("#6ed86e");
	public static final Color DEFAULT_FILL_COLOR = Color.decode("#bcffdd");
	public static final double DEFAULT_RADIUS = 16.0;
	
	private ArrayList<ModeleBase> _bases;
	private Color _outlineColor = DEFAULT_OUTLINE_COLOR;
	private Color _fillColor = DEFAULT_FILL_COLOR;
	private double _radius = DEFAULT_RADIUS;

  public HighlightRegionAnnotation(RNA r, int startIndex, int stopIndex)
  {
	  this(r.getBasesBetween(startIndex, stopIndex));
  }

  public HighlightRegionAnnotation(ArrayList<ModeleBase> b)
  {
	  this(b,DEFAULT_FILL_COLOR,DEFAULT_OUTLINE_COLOR,DEFAULT_RADIUS);
  }
 
  
  public HighlightRegionAnnotation(ArrayList<ModeleBase> b,Color fill, Color outline, double radius)
  {
	  _bases = b;
	  _fillColor = fill;
	  _outlineColor = outline;
	  _radius = radius;
  }
  
	public HighlightRegionAnnotation clone()
	{
		return new HighlightRegionAnnotation(_bases,_fillColor,_outlineColor,_radius);
	}
	
  public int getMinIndex()
  {
	  int min = Integer.MAX_VALUE;
	  for (ModeleBase mb : _bases)
	  {
		  min = Math.min(min, mb.getIndex());
	  }
	  return min;
  }

  public int getMaxIndex()
  {
	  int max = Integer.MIN_VALUE;
	  for (ModeleBase mb : _bases)
	  {
		  max = Math.max(max, mb.getIndex());
	  }
	  return max;
  }

  
  public void setOutlineColor(Color c)
  {
	  _outlineColor = c;  
  }

  public ArrayList<ModeleBase> getBases()
  {
	  return _bases;
  }

  public void setBases(ArrayList<ModeleBase> b)
  {
	  _bases = b;
  }
  
  public void setFillColor(Color c)
  {
	  _fillColor = c;  
  }

  public Color getFillColor()
  {
	  return _fillColor;  
  }
 
  public Color getOutlineColor()
  {
	  return _outlineColor;  
  }

  public double getRadius()
  {
	  return _radius;  
  }

  public void setRadius(double v)
  {
	  _radius = v;  
  }
  
  public static final int NUM_STEPS_ROUNDED_CORNERS = 10;
	
	public GeneralPath getShape(Point2D.Double[] realCoords,Point2D.Double[] realCenters, double scaleFactor)
	{
		GeneralPath p = new GeneralPath();
		LinkedList<Point2D.Double> pointList = new LinkedList<Point2D.Double>();
		boolean isDirect = true;
		if (getBases().size()>=2)
		{
			int j = getBases().get(0).getIndex();
			Point2D.Double p1 = realCoords[j];
			Point2D.Double p2 = realCoords[j+1];
			Point2D.Double p3 = realCenters[j];
			isDirect = RNA.testDirectionality(p1, p3, p2);
		}
		if (getBases().size()>0)
		{
			int j = getBases().get(0).getIndex();
			Point2D.Double point = realCoords[j];
			Point2D.Double centerBck = realCenters[j];
			double dist = point.distance(centerBck);
			Point2D.Double vn = new Point2D.Double((centerBck.x-point.x)/dist,(centerBck.y-point.y)/dist);
			for(int k=1;k<=NUM_STEPS_ROUNDED_CORNERS;k++)
			{
				double angle = k*Math.PI/((double)NUM_STEPS_ROUNDED_CORNERS+1);
				if (isDirect)
					angle += Math.PI;
				Point2D.Double nvn = new Point2D.Double(Math.cos(angle)*vn.x+Math.sin(angle)*vn.y,-Math.sin(angle)*vn.x+Math.cos(angle)*vn.y);
				Point2D.Double interForward = new Point2D.Double(point.x + scaleFactor *getRadius()*nvn.x, 
						point.y + scaleFactor *getRadius()*nvn.y);
				pointList.addLast(interForward);
			}
			
		}
		for (int i = 0;i<getBases().size();i++)
		{ 
			int j1 = getBases().get(i).getIndex();
			if ((j1>0)&&(j1<realCoords.length-1))
			{
				int j0 = j1-1;
				int j2 = j1+1;

				Point2D.Double p0 = realCoords[j0];
				Point2D.Double p1 = realCoords[j1];			
				Point2D.Double p2 = realCoords[j2];

				double dist1 = p2.distance(p1);
				Point2D.Double v1 = new Point2D.Double((p2.x-p1.x)/dist1,(p2.y-p1.y)/dist1);
				Point2D.Double vn1 = new Point2D.Double(v1.y,-v1.x);
				double dist2 = p1.distance(p0);
				Point2D.Double v2 = new Point2D.Double((p1.x-p0.x)/dist2,(p1.y-p0.y)/dist2);
				Point2D.Double vn2 = new Point2D.Double(v2.y,-v2.x);
				double h = (new Point2D.Double(vn2.x-vn1.x,vn2.y-vn1.y).distance(new Point2D.Double(0,0))/2.0);
				Point2D.Double vn = new Point2D.Double((vn1.x+vn2.x)/2.0,(vn1.y+vn2.y)/2.0);
				double D = vn.distance(new Point2D.Double(0.0,0.0));
				vn.x/= D;
				vn.y/= D;
				double nnorm = (D+h*h/D);
				
				double nnormF = nnorm;
				double nnormB = nnorm;
				
				
				Point2D.Double interForward = new Point2D.Double(p1.x + nnormF*scaleFactor *getRadius()*vn.x, 
						p1.y + nnormF*scaleFactor *getRadius()*vn.y);
				Point2D.Double interBackward = new Point2D.Double(p1.x - nnormB*scaleFactor *getRadius()*vn.x, 
						p1.y - nnormB*scaleFactor *getRadius()*vn.y);


				if (pointList.size()>0)
				{
					Point2D.Double prev1 = pointList.getLast();			
					Point2D.Double prev2 = pointList.getFirst();

					if ((interForward.distance(prev1)+interBackward.distance(prev2))<(interForward.distance(prev2)+interBackward.distance(prev1)))
					{
						pointList.addLast(interForward);
						pointList.addFirst(interBackward);
					}
					else
					{
						pointList.addFirst(interForward);
						pointList.addLast(interBackward);
					}
				}
				else
				{
					pointList.addLast(interForward);
					pointList.addFirst(interBackward);
				}
			}
		}
		if (getBases().size()>=2)
		{
			int j = getBases().get(getBases().size()-1).getIndex();
			Point2D.Double p1 = realCoords[j];
			Point2D.Double p2 = realCoords[j-1];
			Point2D.Double p3 = realCenters[j];
			isDirect = RNA.testDirectionality(p1, p2, p3);
		}

		if (getBases().size()>0)
		{
			int j = getBases().get(getBases().size()-1).getIndex();
			Point2D.Double point = realCoords[j];
			Point2D.Double centerBck = realCenters[j];
			double dist = point.distance(centerBck);
			Point2D.Double vn = new Point2D.Double((centerBck.x-point.x)/dist,(centerBck.y-point.y)/dist);
			for(int k=1;k<=NUM_STEPS_ROUNDED_CORNERS;k++)
			{
				double angle = k*Math.PI/((double)NUM_STEPS_ROUNDED_CORNERS+1);
				if (!isDirect)
					angle += Math.PI;
				Point2D.Double nvn = new Point2D.Double(Math.cos(angle)*vn.x+Math.sin(angle)*vn.y,-Math.sin(angle)*vn.x+Math.cos(angle)*vn.y);
				Point2D.Double interForward = new Point2D.Double(point.x + scaleFactor *getRadius()*nvn.x, 
						point.y + scaleFactor *getRadius()*nvn.y);
				pointList.addLast(interForward);
			}
			
		}
		

		if (pointList.size()>0)
		{
			Point2D.Double point = pointList.get(0);
			p.moveTo((float)point.x, (float)point.y);

			for (int i=1;i<pointList.size();i++)
			{
				point = pointList.get(i);
				p.lineTo((float)point.x, (float)point.y);				
			}
			p.closePath();
		}
		return p;
	}

	
	public static HighlightRegionAnnotation parseHighlightRegionAnnotation(String txt, VARNAPanel vp)
	{
		try
		{
		String[] parts = txt.split(":");
		String[] coords = parts[0].split("-");
		int from = Integer.parseInt(coords[0]);
		int to = Integer.parseInt(coords[1]);
		int  i = vp.getRNA().getIndexFromBaseNumber(from);
		int  j = vp.getRNA().getIndexFromBaseNumber(to);
		Color fill = HighlightRegionAnnotation.DEFAULT_FILL_COLOR;
		Color outline = HighlightRegionAnnotation.DEFAULT_OUTLINE_COLOR;
		double radius = HighlightRegionAnnotation.DEFAULT_RADIUS;
		ArrayList<ModeleBase> bases = vp.getRNA().getBasesBetween(i, j);
		try
		{
			String[] options = parts[1].split(",");
			for (int k = 0; k < options.length; k++) 
			{
				//System.out.println(options[k]);
				try
				{
					String[] data = options[k].split("=");
					String lhs = data[0].toLowerCase();
					String rhs = data[1];
					if (lhs.equals("fill"))
					{
						fill = VARNAConfigLoader.getSafeColor(rhs, fill);
					}
					else if (lhs.equals("outline"))
					{
						outline = VARNAConfigLoader.getSafeColor(rhs, outline);
					}
					else if (lhs.equals("radius"))
					{
						radius = Double.parseDouble(rhs);
					}	
				}
				catch(Exception e)
				{
				}				
			}
		}
		catch(Exception e)
		{
		}
		return new HighlightRegionAnnotation(bases,fill,outline,radius);
		}
	catch(Exception e)
	{
	}
	return null;
	}
	
	public String toString()
	{
		//String result = "HighlightRegionAnnotation[";
		//result += "fill:"+_fillColor.toString();
		//result += ",outline:"+_outlineColor.toString();
		//result += ",radius:"+_radius;
		//return result+"]";
		String result = "Highlighted region "+getMinIndex()+"-"+getMaxIndex();
		return result;
	}
	
}
