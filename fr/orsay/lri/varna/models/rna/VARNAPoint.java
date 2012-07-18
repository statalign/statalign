package fr.orsay.lri.varna.models.rna;

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;
import java.io.Serializable;

public class VARNAPoint implements Serializable {

	private static final long serialVersionUID = 8815373295131046029L;
	
	public double x = 0.0;
	public double y = 0.0;
	public VARNAPoint()
	{ this(0.0,0.0); }
	public VARNAPoint(double px, double py)
    {
    	x = px; y = py;
    }
    public VARNAPoint(Point2D.Double p)
    {
    	this(p.x,p.y);
    }
 
    public double getX()
    {
    	return x;
    }
    
    public double getY()
    {
    	return y;
    }
    
    public Point2D.Double toPoint2D()
    {
    	return new Point2D.Double(x,y);
    }
    
    public String toString()
    {
    	return "("+x+","+y+")" ;
    }
}
