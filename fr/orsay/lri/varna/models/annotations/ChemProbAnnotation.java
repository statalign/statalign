package fr.orsay.lri.varna.models.annotations;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.Serializable;

import fr.orsay.lri.varna.models.rna.ModeleBase;


public class ChemProbAnnotation implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 5833315460145031242L;


	public enum ChemProbAnnotationType 
	{
		TRIANGLE_TYPE,
		ARROW_TYPE,
		PIN_TYPE,
		DOT_TYPE;
		public String toString()
		{
			switch(this)
			{
			case TRIANGLE_TYPE:
				return "Triangle";
			case ARROW_TYPE:
				return "Arrow";
			case DOT_TYPE:
				return "Dot";
			case PIN_TYPE:
				return "Pin";
			default:
				return "Unknown";
			}
		}
	};
	
	public static double DEFAULT_INTENSITY = 1.0;
	public static ChemProbAnnotationType DEFAULT_TYPE = ChemProbAnnotationType.ARROW_TYPE;
	public static Color DEFAULT_COLOR = Color.blue.darker();

	private ModeleBase _mbfst;
	private ModeleBase _mbsnd;
	private Color _color;
	private double _intensity;
	private ChemProbAnnotationType _type;
	private boolean _out;
	
	public ChemProbAnnotation(ModeleBase mbfst, ModeleBase mbsnd, String styleDesc) {
		this(mbfst,mbsnd);
		applyStyle(styleDesc);
	}

	public ChemProbAnnotation(ModeleBase mbfst, ModeleBase mbsnd) {
		this(mbfst,mbsnd,ChemProbAnnotation.DEFAULT_TYPE,ChemProbAnnotation.DEFAULT_INTENSITY);
	}

	public ChemProbAnnotation(ModeleBase mbfst, ModeleBase mbsnd, double intensity) {
		this(mbfst,mbsnd,ChemProbAnnotation.DEFAULT_TYPE,intensity);
	}

	public ChemProbAnnotation(ModeleBase mbfst, ModeleBase mbsnd, ChemProbAnnotationType type) {
		this(mbfst,mbsnd,type,ChemProbAnnotation.DEFAULT_INTENSITY);
	}
	
	public ChemProbAnnotation(ModeleBase mbfst, ModeleBase mbsnd, ChemProbAnnotationType type, double intensity) {
		if (mbfst.getIndex()>mbsnd.getIndex())
		{
			ModeleBase tmp = mbsnd;
			mbsnd = mbfst;
			mbfst = tmp;
		}
		_mbfst = mbfst;
		_mbsnd = mbsnd;
		_type = type;
		_intensity = intensity;
		_color = DEFAULT_COLOR;
		_out = true;
	}

	public boolean isOut()
	{
		return _out; 
	}
	
	public void setOut(boolean b)
	{
		_out = b;
	}

	public Color getColor()
	{
		return _color; 
	}

	public double getIntensity()
	{
		return _intensity; 
	}
	
	public ChemProbAnnotationType getType()
	{
		return _type; 
	}
	
	public void setColor(Color c){
		_color = c;
	}

	public void setIntensity(double d){
		_intensity = d;
	}
	
	public Point2D.Double getAnchorPosition()
	{
		Point2D.Double result = new Point2D.Double(
				(_mbfst.getCoords().x+_mbsnd.getCoords().x)/2.0,
				(_mbfst.getCoords().y+_mbsnd.getCoords().y)/2.0);
		return result;
	}
	
	public Point2D.Double getDirVector()
	{
		Point2D.Double norm = getNormalVector();
		Point2D.Double result = new Point2D.Double(-norm.y,norm.x);
		Point2D.Double anchor = getAnchorPosition();
		Point2D.Double center = new Point2D.Double(
				(_mbfst.getCenter().x+_mbsnd.getCenter().x)/2.0,
				(_mbfst.getCenter().y+_mbsnd.getCenter().y)/2.0);
		Point2D.Double vradius = new Point2D.Double(
				(center.x-anchor.x)/2.0,
				(center.y-anchor.y)/2.0);

		if (_out)
		{
		  if (result.x*vradius.x+result.y*vradius.y>0)
		  {
			return new Point2D.Double(-result.x,-result.y);
		  }
		}
		else
		{
			  if (result.x*vradius.x+result.y*vradius.y<0)
			  {
				return new Point2D.Double(-result.x,-result.y);
			  }			
		}
		return result;		
	}
	public Point2D.Double getNormalVector()
	{
		Point2D.Double tmp = new Point2D.Double(
				(_mbsnd.getCoords().x-_mbfst.getCoords().x)/2.0,
				(_mbsnd.getCoords().y-_mbfst.getCoords().y)/2.0);
		double norm = tmp.distance(0, 0);
		Point2D.Double result = new Point2D.Double(tmp.x/norm,tmp.y/norm);
		return result;				
	}
	
	
	public void applyStyle(String styleDesc)
	{
		String[] chemProbs = styleDesc.split(",");
		for (int i = 0; i < chemProbs.length; i++) {
			String thisStyle = chemProbs[i];
			String[] data = thisStyle.split("=");
			if (data.length==2)
			{
				String name = data[0];
				String value = data[1];
				if (name.toLowerCase().equals("color"))
				{
					Color c = Color.decode(value);
					if (c==null)
					{ c = _color; }
					setColor(c);
				}
				else if (name.toLowerCase().equals("intensity"))
				{
					_intensity = Double.parseDouble(value);
				}
				else if (name.toLowerCase().equals("dir"))
				{
					_out = value.toLowerCase().equals("out");
				}
				else if (name.toLowerCase().equals("glyph"))
				{
					if (value.toLowerCase().equals("arrow"))
					{_type = ChemProbAnnotationType.ARROW_TYPE;}
					else if (value.toLowerCase().equals("triangle"))
					{_type = ChemProbAnnotationType.TRIANGLE_TYPE;}
					else if (value.toLowerCase().equals("pin"))
					{_type = ChemProbAnnotationType.PIN_TYPE;}
					else if (value.toLowerCase().equals("dot"))
					{_type = ChemProbAnnotationType.DOT_TYPE;}
				}
			}
		}
	}
	
	public void setType(ChemProbAnnotationType s)
	{
		_type = s;
	}

	public ChemProbAnnotation clone()
	{
		ChemProbAnnotation result = new ChemProbAnnotation(this._mbfst,this._mbsnd);
		result._intensity = _intensity;
		result._type = _type;
		result._color= _color;
		result._out = _out;
		return result;
	}
	
	public String toString()
	{
		return "Chem. prob. "+this._type+" Base#"+this._mbfst.getBaseNumber()+"-"+this._mbsnd.getBaseNumber();
	}
	
}
