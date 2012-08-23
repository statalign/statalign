/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Université Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.models.rna;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Point;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Reader;
import java.io.Serializable;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.InputSource;

import statalign.postprocess.gui.ColorGradient;
import statalign.postprocess.gui.StructureGUI;
import statalign.postprocess.plugins.Structure;

import fr.orsay.lri.varna.exceptions.ExceptionExportFailed;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionInvalidRNATemplate;
import fr.orsay.lri.varna.exceptions.ExceptionLoadingFailed;
import fr.orsay.lri.varna.exceptions.ExceptionNAViewAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionPermissionDenied;
import fr.orsay.lri.varna.exceptions.ExceptionUnmatchedClosingParentheses;
import fr.orsay.lri.varna.exceptions.ExceptionWritingForbidden;
import fr.orsay.lri.varna.interfaces.InterfaceVARNAListener;
import fr.orsay.lri.varna.interfaces.InterfaceVARNAObservable;
import fr.orsay.lri.varna.models.CubicBezierCurve;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.annotations.ChemProbAnnotation;
import fr.orsay.lri.varna.models.annotations.HighlightRegionAnnotation;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.export.PSExport;
import fr.orsay.lri.varna.models.export.SVGExport;
import fr.orsay.lri.varna.models.export.SecStrDrawingProducer;
import fr.orsay.lri.varna.models.export.XFIGExport;
import fr.orsay.lri.varna.models.naView.NAView;
import fr.orsay.lri.varna.models.templates.RNANodeValue2TemplateDistance;
import fr.orsay.lri.varna.models.templates.RNANodeValueTemplate;
import fr.orsay.lri.varna.models.templates.RNANodeValueTemplateBasePair;
import fr.orsay.lri.varna.models.templates.RNATemplate;
import fr.orsay.lri.varna.models.templates.RNATemplateAlign;
import fr.orsay.lri.varna.models.templates.RNATemplateDrawingAlgorithmException;
import fr.orsay.lri.varna.models.templates.RNATemplateMapping;
import fr.orsay.lri.varna.models.templates.RNATemplate.EdgeEndPointPosition;
import fr.orsay.lri.varna.models.templates.RNATemplate.In1Is;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateHelix;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateUnpairedSequence;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement.EdgeEndPoint;
import fr.orsay.lri.varna.models.treealign.RNANodeValue2;
import fr.orsay.lri.varna.models.treealign.RNATree2;
import fr.orsay.lri.varna.models.treealign.RNATree2Exception;
import fr.orsay.lri.varna.models.treealign.Tree;
import fr.orsay.lri.varna.models.treealign.TreeAlign;
import fr.orsay.lri.varna.models.treealign.TreeAlignException;
import fr.orsay.lri.varna.models.treealign.TreeAlignResult;
import fr.orsay.lri.varna.views.VueUI;

import java.lang.Math;


/**
 * The RNA model which contain the base list and the draw algorithm mode
 * 
 * @author darty
 * 
 */
public class RNA extends InterfaceVARNAObservable implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 7541274455751497303L;


	/**
	 * Selects the "Feynman diagram" drawing algorithm that places the bases on
	 * a circle and draws the base-pairings as chords of the circle graph.
	 */
	
	public static final int DRAW_MODE_CIRCULAR = 1;
	/**
	 * Selects the "tree drawing" algorithm. Draws each loop on a circle whose
	 * radius depends on the number of bases involved in the loop. As some
	 * helices can be overlapping in the result, basic interaction is provided
	 * so that the user can "disentangle" the drawing by spinning the helices
	 * around the axis defined by their multiloop (bulge or internal loop)
	 * origin. This is roughly the initial placement strategy of RNAViz.
	 * 
	 * @see <a href="http://rnaviz.sourceforge.net/">RNAViz</a>
	 */
	public static final int DRAW_MODE_RADIATE = 2;
	
	/**
	 * Selects the NAView algorithm.
	 */
	public static final int DRAW_MODE_NAVIEW = 3;
	/**
	 * Selects the linear algorithm.
	 */
	public static final int DRAW_MODE_LINEAR = 4;

	public static final int DRAW_MODE_VARNA_VIEW = 5;
	
	/**
	 * Selects the RNAView algorithm.
	 */
	public static final int DRAW_MODE_MOTIFVIEW = 6;
	
	public static final int DRAW_MODE_TEMPLATE = 7;

	public static final int DEFAULT_DRAW_MODE = DRAW_MODE_RADIATE;
	
	private Double _spaceBetweenBases = 1.0;
	/**
	 * The draw algorithm mode
	 */
	private int _drawMode = DRAW_MODE_RADIATE;

	public int BASE_RADIUS = 10;
	public static final double LOOP_DISTANCE = 40.0; // distance between base pairs in an helix
	public static final double BASE_PAIR_DISTANCE = 65.0; // distance between the two bases of a pair
	public static final double MULTILOOP_DISTANCE = 35.0;
	public static final double VIRTUAL_LOOP_RADIUS= 40.0;

	public double CHEM_PROB_DIST = 14;
	public double CHEM_PROB_BASE_LENGTH = 30;
	public double CHEM_PROB_ARROW_HEIGHT = 10;
	public double CHEM_PROB_ARROW_WIDTH = 5;
	public double CHEM_PROB_TRIANGLE_WIDTH = 2.5;
	public double CHEM_PROB_PIN_SEMIDIAG = 6;
	public double CHEM_PROB_DOT_RADIUS = 6.;
	
	public GeneralPath _debugShape = null;
	
	
	private boolean _drawn = false;
	public boolean probMode = true;
	
	private ColorGradient cg = new ColorGradient(Color.WHITE, Color.BLUE);
	
	public double _bpHeightIncrement = VARNAConfig.DEFAULT_BP_INCREMENT;

	/**
	 * the base list
	 */
	protected ArrayList<ModeleBase> _listeBases;
	
	/**
	 * the strand list
	 */
	StructureTemp _listStrands = new StructureTemp();
	
	/**
	 * Additional bonds and info can be specified here.
	 */
	private ArrayList<ModeleBP> _structureAux = new ArrayList<ModeleBP>();

	transient private ArrayList<InterfaceVARNAListener> _listeVARNAListener = new ArrayList<InterfaceVARNAListener>();

	static ArrayList<String> _normalBases = new ArrayList<String>();
	{
		_normalBases.add("a");
		_normalBases.add("c");
		_normalBases.add("g");
		_normalBases.add("u");
		_normalBases.add("t");
	}
	
	private ArrayList<TextAnnotation> _listeAnnotations = new ArrayList<TextAnnotation>();
	private ArrayList<HighlightRegionAnnotation> _listeRegionHighlights = new ArrayList<HighlightRegionAnnotation>();

	
	private String _name = "";

	public RNA() {
		this("");
	}
	public RNA(String n) {
		_name = n;
		_listeBases = new ArrayList<ModeleBase>();
		_drawn = false;
		init();
	}

	public String toString()
	{
		if (_name.equals(""))
		{
			return getStructDBN();
		}
		else
		{
			return _name; 
		}
	}
	
	public RNA(RNA r) {
		_spaceBetweenBases = r._spaceBetweenBases;
		_drawMode = r._drawMode;
		_listeBases.addAll(r._listeBases);
		_listeVARNAListener = (ArrayList<InterfaceVARNAListener>) r._listeVARNAListener;
		_drawn = r._drawn;
		init();
	}

	public void init() {
	}

	public void saveRNADBN(String path, String title)
			throws ExceptionWritingForbidden {
		try {
			FileWriter out = new FileWriter(path);
			if (!title.equals("")) {
				out.write("> " + title + "\n");
			}
			out.write(getListeBasesToString());
			out.write('\n');
			String str = "";
			for (int i = 0; i < _listeBases.size(); i++) {
				if (_listeBases.get(i).getElementStructure() == -1) {
					str += '.';
				} else {
					if (_listeBases.get(i).getElementStructure() > i) {
						str += '(';
					} else {
						str += ')';
					}
				}
			}
			out.write(str);
			out.write('\n');
			out.close();
		} catch (IOException e) {
			throw new ExceptionWritingForbidden(e.getMessage());
		}
	}

	public Color getBaseInnerColor(int i, VARNAConfig conf) {
		/*Color result = _listeBases.get(i).getStyleBase()
				.get_base_inner_color();
		String res = _listeBases.get(i).getContent();
		if (conf._drawColorMap)
		{
			result = conf._cm.getColorForValue(_listeBases.get(i).getValue());			
		}
		else if ((conf._colorDashBases && (res.contains("-")))) {
			result = conf._dashBasesColor;
		}
		else if ((conf._colorSpecialBases && !_normalBases.contains(res.toLowerCase()))) {
			result = conf._specialBasesColor;
		} */
		
		Color result = _listeBases.get(i).getStyleBase()
		.get_base_inner_color();
		String res = _listeBases.get(i).getContent();
		
		if(!probMode) {
			if(res.toLowerCase().equals("a")) {result = Color.RED;}
			
			else if(res.toLowerCase().equals("c")) {result = Color.BLUE;}
			else if(res.toLowerCase().equals("g")) {result = Color.YELLOW;}
			else if(res.toLowerCase().equals("u") || res.toLowerCase().equals("t")) {result = Color.GREEN;}
			
			else if ((conf._colorSpecialBases && !_normalBases.contains(res.toLowerCase()))) {
				result = conf._specialBasesColor;
			} 
		}
		
		else {
			if(Structure.currentDotBracketStructure.charAt(i) == '(') {
				int j = _listeBases.get(i).getStyleBP().getIndex3();
				result = cg.getColor(Structure.probMatrix[i][j]);
			}
			
			else if(Structure.currentDotBracketStructure.charAt(i) == ')') {
				int j = _listeBases.get(i).getStyleBP().getIndex5();
				result = cg.getColor(Structure.probMatrix[i][j]);
			}
			
			else if(Structure.currentDotBracketStructure.charAt(i) == '.') {
				result = cg.getColor(Structure.singleMatrix[i]);
			}
		}
		
		return result;
	}

	public Color getBaseOuterColor(int i, VARNAConfig conf) {
		Color result = _listeBases.get(i).getStyleBase()
				.get_base_outline_color();
		return result;
	}

	public Color getBaseNameColor(int i, VARNAConfig conf) {
		Color result = _listeBases.get(i).getStyleBase().get_base_name_color();
		return result;
	}

	public Color getBasePairColor(ModeleBP bp, VARNAConfig conf) {
		Color bondColor = conf._bondColor;
		if (conf._useBaseColorsForBPs) {
			bondColor = _listeBases.get(bp.getPartner5().getIndex())
					.getStyleBase().get_base_inner_color();
		}
		if (bp != null) {
			bondColor = bp.getStyle().getColor(bondColor);
		}
		return bondColor;
	}

	public double getBasePairThickness(ModeleBP bp, VARNAConfig conf) {
		double thickness = bp.getStyle().getThickness(conf._bpThickness);
		return thickness;
	}

	private void drawSymbol(SecStrDrawingProducer out, double posx,
			double posy, double normx, double normy, double radius,
			boolean isCIS, ModeleBP.Edge e, double thickness) {
		Color bck = out.getCurrentColor();
		switch (e) {
		case WATSON_CRICK:
			if (isCIS) {
				out.fillCircle(posx, posy, (radius / 2.0), thickness, bck);
			} else {
				out.drawCircle(posx, posy, (radius / 2.0), thickness);
				out.fillCircle(posx, posy, (radius / 2.0), thickness,
						Color.white);
			}
			break;
		case HOOGSTEEN: {
			double xtab[] = new double[4];
			double ytab[] = new double[4];
			xtab[0] = posx - radius * normx / 2.0 - radius * normy / 2.0;
			ytab[0] = posy - radius * normy / 2.0 + radius * normx / 2.0;
			xtab[1] = posx + radius * normx / 2.0 - radius * normy / 2.0;
			ytab[1] = posy + radius * normy / 2.0 + radius * normx / 2.0;
			xtab[2] = posx + radius * normx / 2.0 + radius * normy / 2.0;
			ytab[2] = posy + radius * normy / 2.0 - radius * normx / 2.0;
			xtab[3] = posx - radius * normx / 2.0 + radius * normy / 2.0;
			ytab[3] = posy - radius * normy / 2.0 - radius * normx / 2.0;
			if (isCIS) {
				out.fillPolygon(xtab, ytab, bck);
			} else {
				out.drawPolygon(xtab, ytab, thickness);
				out.fillPolygon(xtab, ytab,  Color.white);
			}
		}
			break;
		case SUGAR: {
			double ix = radius * normx / 2.0;
			double iy = radius * normy / 2.0;
			double jx = radius * normy / 2.0;
			double jy = -radius * normx / 2.0;
			double xtab[] = new double[3];
			double ytab[] = new double[3];
			xtab[0] = posx - ix + jx;
			ytab[0] = posy - iy + jy;
			xtab[1] = posx + ix + jx;
			ytab[1] = posy + iy + jy;
			xtab[2] = posx - jx;
			ytab[2] = posy - jy;

			if (isCIS) {
				out.fillPolygon(xtab, ytab, bck);
			} else {
				out.drawPolygon(xtab, ytab, thickness);
				out.fillPolygon(xtab, ytab,  Color.white);
			}
		}
			break;
		}
	}

	private void drawBasePair(SecStrDrawingProducer out, Point2D.Double orig,
			Point2D.Double dest, ModeleBP style, VARNAConfig conf) {
		double dx = dest.x - orig.x;
		double dy = dest.y - orig.y;
		double dist = Math.sqrt((dest.x - orig.x) * (dest.x - orig.x)
				+ (dest.y - orig.y) * (dest.y - orig.y));
		dx /= dist;
		dy /= dist;
		double nx = -dy;
		double ny = dx;
		orig = new Point2D.Double(orig.x+BASE_RADIUS*dx,orig.y+BASE_RADIUS*dy);
		dest = new Point2D.Double(dest.x-BASE_RADIUS*dx,dest.y-BASE_RADIUS*dy);
		if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_LW) {
			double thickness = getBasePairThickness(style, conf);
			double radiusCircle = ((BASE_PAIR_DISTANCE - BASE_RADIUS) / 5.0);

			if (style.isCanonical()) {
				if (style.isCanonicalGC()) {
					if ((orig.x != dest.x) || (orig.y != dest.y)) {
						nx *= BASE_RADIUS / 4.0;
						ny *= BASE_RADIUS / 4.0;
						out
								.drawLine((orig.x + nx), (orig.y + ny),
										(dest.x + nx), (dest.y + ny),
										conf._bpThickness);
						out
								.drawLine((orig.x - nx), (orig.y - ny),
										(dest.x - nx), (dest.y - ny),
										conf._bpThickness);
					}
				} else if (!style.isWobbleUG()) {
					double cx = (dest.x + orig.x) / 2.0;
					double cy = (dest.y + orig.y) / 2.0;
					out.drawLine(orig.x, orig.y, dest.x, dest.y,
							conf._bpThickness);
					drawSymbol(out, cx, cy, nx, ny, radiusCircle,
							style.isCIS(), style.getEdgePartner5(), thickness);
				} else {
					out.drawLine(orig.x, orig.y, dest.x, dest.y,
							conf._bpThickness);
				}
			} else {
				ModeleBP.Edge p1 = style.getEdgePartner5();
				ModeleBP.Edge p2 = style.getEdgePartner3();
				double cx = (dest.x + orig.x) / 2.0;
				double cy = (dest.y + orig.y) / 2.0;
				out.drawLine(orig.x, orig.y, dest.x, dest.y, conf._bpThickness);
				if (p1 == p2) {
					drawSymbol(out, cx, cy, nx, ny, radiusCircle,
							style.isCIS(), p1, thickness);
				} else {
					double vdx = (dest.x - orig.x);
					double vdy = (dest.y - orig.y);
					vdx /= 6.0;
					vdy /= 6.0;
					drawSymbol(out, cx + vdx, cy + vdy, nx, ny, radiusCircle,
							style.isCIS(), p2, thickness);
					drawSymbol(out, cx - vdx, cy - vdy, nx, ny, radiusCircle,
							style.isCIS(), p1, thickness);
				}
			}
		} else if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_RNAVIZ) {
			double xcenter = (orig.x + dest.x) / 2.0;
			double ycenter = (orig.y + dest.y) / 2.0;
			out.fillCircle(xcenter, ycenter, 3.0 * conf._bpThickness,
					conf._bpThickness, out.getCurrentColor());
		} else if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_SIMPLE) {
		}
	}

	
	private void drawColorMap(VARNAConfig _conf,SecStrDrawingProducer out)
	{
		double v1 = _conf._cm.getMinValue();
		double v2 = _conf._cm.getMaxValue();
		int x,y;
		double xSpaceAvail = 0;
		double ySpaceAvail = 0;
		double thickness = 1.0;
		/*ySpaceAvail = Math.min((getHeight()-rnabbox.height*scaleFactor-getTitleHeight())/2.0,scaleFactor*(_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE));
		if ((int)ySpaceAvail==0)
		{
			xSpaceAvail = Math.min((getWidth()-rnabbox.width*scaleFactor)/2,scaleFactor*(_conf._colorMapWidth)+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH);			
		}*/
		Rectangle2D.Double currentBBox = out.getBoundingBox(); 
		
		double xBase =  (currentBBox.getMaxX() -_conf._colorMapWidth-_conf._colorMapXOffset);
		//double yBase =  (minY - _conf._colorMapHeight + _conf._colorMapYOffset);
		double yBase =  (currentBBox.getMinY()-_conf._colorMapHeight-VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE);
		
		for (int i=0;i<_conf._colorMapWidth;i++)
		{
			double ratio = (((double)i)/((double)_conf._colorMapWidth-1));
			double val = v1+(v2-v1)*ratio;
			Color c = _conf._cm.getColorForValue(val);
			x = (int) (xBase + i);
			y = (int) yBase;
			out.fillRectangle(x, y, VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH, _conf._colorMapHeight,c);	
		}
		out.setColor(VARNAConfig.DEFAULT_COLOR_MAP_OUTLINE);
		out.drawRectangle(xBase,yBase, (double)_conf._colorMapWidth+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH-1, _conf._colorMapHeight,thickness);
		
		out.setColor(VARNAConfig.DEFAULT_COLOR_MAP_FONT_COLOR);
		out.setFont(out.getCurrentFont(),VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.5);
		out.drawText(xBase,
				yBase+_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7,
				""+_conf._cm.getMinValue());
		out.drawText(xBase+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH+_conf._colorMapWidth,
				yBase+_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7,
				""+_conf._cm.getMaxValue());
			out.drawText(xBase+(VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH+_conf._colorMapWidth)/2.0,
				yBase-(VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7),
				_conf._colorMapCaption);
		
	}

	
	private void renderRegionHighlights(SecStrDrawingProducer out, Point2D.Double[] realCoords,Point2D.Double[] realCenters)
	{
		for (HighlightRegionAnnotation r:_listeRegionHighlights)
		{
			GeneralPath s = r.getShape(realCoords,realCenters,1.0);
			out.setColor(r.getFillColor());
			out.fillPolygon(s, r.getFillColor());
			out.setColor(r.getOutlineColor());
			out.drawPolygon(s, 1l);
		}

	}
	
	private void saveRNA(String path, VARNAConfig conf, double scale,
			SecStrDrawingProducer out) throws ExceptionWritingForbidden {
		out.setScale(scale);
		// Computing bounding boxes
		double EPSMargin = 40;
		double minX = Double.MAX_VALUE;
		double maxX = Double.MIN_VALUE;
		double minY = Double.MAX_VALUE;
		double maxY = Double.MIN_VALUE;

		double x0, y0, x1, y1, xc, yc, xp, yp, dx, dy, norm;

		for (int i = 0; i < _listeBases.size(); i++) {
			minX = Math.min(minX, (_listeBases.get(i).getCoords().getX()
					- BASE_RADIUS - EPSMargin));
			minY = Math.min(minY, -(_listeBases.get(i).getCoords().getY()
					- BASE_RADIUS - EPSMargin));
			maxX = Math.max(maxX, (_listeBases.get(i).getCoords().getX()
					+ BASE_RADIUS + EPSMargin));
			maxY = Math.max(maxY, -(_listeBases.get(i).getCoords().getY()
					+ BASE_RADIUS + EPSMargin));
		}
		
		// Rescaling everything
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()]; 
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()]; 
		for (int i = 0; i < _listeBases.size(); i++) {
			xp = (_listeBases.get(i).getCoords().getX() - minX);
			yp = -(_listeBases.get(i).getCoords().getY() - minY);
			coords[i] = new Point2D.Double(xp,yp);

			Point2D.Double centerBck = getCenter(i);
			if (get_drawMode() == RNA.DRAW_MODE_NAVIEW
					|| get_drawMode() == RNA.DRAW_MODE_RADIATE) 
			{
				if ((_listeBases.get(i).getElementStructure() != -1)
						&& i < _listeBases.size() - 1
						&& i > 1)
				{
				  ModeleBase b1 = get_listeBases().get(i - 1);
				  ModeleBase b2 = get_listeBases().get(i + 1);
				  int j1 = b1.getElementStructure();
				  int j2 = b2.getElementStructure();
				  if ((j1==-1)^(j2==-1))
				  {
				  // alors la position du nombre associé doit etre décalé
					Point2D.Double a1 = b1.getCoords();
					Point2D.Double a2 = b2.getCoords();
					Point2D.Double c1 = b1.getCenter();
					Point2D.Double c2 = b2.getCenter();
					
					centerBck.x = _listeBases.get(i).getCoords().x + (c1.x-a1.x)/c1.distance(a1)+(c2.x-a2.x)/c2.distance(a2);
					centerBck.y = _listeBases.get(i).getCoords().y + (c1.y-a1.y)/c1.distance(a1)+(c2.y-a2.y)/c2.distance(a2);
				  }
				}
			}
			xc = (centerBck.getX() - minX);
			yc = -(centerBck.getY() - minY);
			centers[i] = new Point2D.Double(xc,yc);
		}


		// Drawing background
		if (conf._drawBackground)
		  out.setBackgroundColor(conf._backgroundColor);
		
		// Drawing region highlights
		renderRegionHighlights(out,coords,centers);
		
		// Drawing backbone
		out.setColor(conf._backboneColor);
		for (int i = 1; i < _listeBases.size(); i++) {
			Point2D.Double p1 = coords[i-1]; 
			Point2D.Double p2 = coords[i]; 
			x0 = p1.x;
			y0 = p1.y;
			x1 = p2.x;
			y1 = p2.y;
			Point2D.Double vn = new Point2D.Double();
			double dist = p1.distance(p2);
			boolean discontinuous= (getBaseNumber(i)-getBaseNumber(i-1)!=1);
			int a = _listeBases.get(i-1).getElementStructure();
			int b = _listeBases.get(i).getElementStructure();
			boolean consecutive = (a==i)&&(b==i-1);

			if (dist>0){
				vn.x = (x1-x0)/dist;
				vn.y = (y1-y0)/dist;
				if (consecutive && (getDrawMode()!=RNA.DRAW_MODE_LINEAR) && (getDrawMode()!=RNA.DRAW_MODE_CIRCULAR))
				{
					int dir = 0;
					if (i+1<coords.length)
					{dir = (testDirectionality(i-1,i,i+1)?-1:1);}
					else if (i-2>=0)
					{dir = (testDirectionality(i-2,i-1,i)?-1:1);}
					Point2D.Double centerSeg = new Point2D.Double((p1.x+p2.x)/2.0,(p1.y+p2.y)/2.0);
					double centerDist = RNA.VIRTUAL_LOOP_RADIUS*scale; 
					Point2D.Double centerLoop = new Point2D.Double(centerSeg.x+centerDist*dir*vn.y,centerSeg.y-centerDist*dir*vn.x);
					out.drawLine(centerLoop.x-5, centerLoop.y,centerLoop.x+5, centerLoop.y, 2.0);
					out.drawLine(centerLoop.x, centerLoop.y-5,centerLoop.x, centerLoop.y+5, 2.0);
					double radius = centerLoop.distance(p1);
					double a1 = 360.*(Math.atan2(dir*(p1.x-centerLoop.x),dir*(p1.y-centerLoop.y)))/(2.*Math.PI);						
					double a2 = 360.*(Math.atan2(dir*(centerLoop.x-p2.x),dir*(centerLoop.y-p2.y)))/(2.*Math.PI);
					if (a1<0){a1+=360.;}
					if (a2<0){a2+=360.;}
					double angle = a2-a1;
					if (dir*(a2-a1)<0.) angle += dir*360.;
					out.drawArc(centerLoop, 2.*radius, 2.*radius, a1, angle);
				}
				else
				{
				out.drawLine((x0+BASE_RADIUS*vn.x),
						(y0+BASE_RADIUS*vn.y),
						(x1-BASE_RADIUS*vn.x),
						(y1-BASE_RADIUS*vn.y),1.0);
				}	
			}
		}

		// Drawing bonds
		for (int i = 0; i < _listeBases.size(); i++) {
			if (_listeBases.get(i).getElementStructure() > i) {
				ModeleBP style = _listeBases.get(i).getStyleBP();
				if (style.isCanonical() || conf._drawnNonCanonicalBP) {
					Color bpcol = getBasePairColor(style, conf);
					out.setColor(bpcol);

					int j = _listeBases.get(i).getElementStructure();
					x0 = coords[i].x;
					y0 = coords[i].y;
					x1 = coords[j].x;
					y1 = coords[j].y;
					dx = x1 - x0;
					dy = y1 - y0;
					norm = Math.sqrt(dx * dx + dy * dy);
					dx /= norm;
					dy /= norm;

					if (_drawMode == DRAW_MODE_CIRCULAR
							|| _drawMode == DRAW_MODE_RADIATE
							|| _drawMode == DRAW_MODE_NAVIEW) {
						drawBasePair(out, new Point2D.Double(x0, y0),
								new Point2D.Double(x1, y1), style, conf);
					} else if (_drawMode == DRAW_MODE_LINEAR) {
						double coef;
						double distance;
						if (j - i == 1)
							coef = _bpHeightIncrement * 2;
						else
							coef = _bpHeightIncrement * 1;
						distance = (int) Math.round(x1 - x0);
						out.drawArc(new Point2D.Double(x0, y0), distance,
								distance * coef, 0, 180);
					}
				}
			}
		}

		// Drawing additional bonds
		if (conf._drawnNonPlanarBP) {
			for (int i = 0; i < _structureAux.size(); i++) {
				ModeleBP bp = _structureAux.get(i);
				out.setColor(getBasePairColor(bp, conf));
				
				int a = bp.getPartner5().getIndex();
				int b = bp.getPartner3().getIndex();

				if (bp.isCanonical() || conf._drawnNonCanonicalBP) {
					x0 = coords[a].x;
					y0 = coords[a].y;
					x1 = coords[b].x;
					y1 = coords[b].y;
					dx = x1 - x0;
					dy = y1 - y0;
					norm = Math.sqrt(dx * dx + dy * dy);
					dx /= norm;
					dy /= norm;
					if ((_drawMode == DRAW_MODE_CIRCULAR)
							|| (_drawMode == DRAW_MODE_RADIATE)
							|| _drawMode == DRAW_MODE_NAVIEW) {
						drawBasePair(out, new Point2D.Double(x0, y0),
								new Point2D.Double(x1, y1), bp, conf);
					} else if (_drawMode == DRAW_MODE_LINEAR) {
						double coef;
						double distance;
						if (b - a == 1)
							coef = _bpHeightIncrement * 2;
						else
							coef = _bpHeightIncrement * 1;
						distance = (int) Math.round(x1 - x0);
						out.drawArc(new Point2D.Double(x0, y0), distance,
								distance * coef, 0, 180);
					}
				}
			}
		}

		// Drawing Bases
		double baseFontSize = (1.5 * BASE_RADIUS);
		out.setFont(PSExport.FONT_HELVETICA_BOLD, baseFontSize);
		if (conf._comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				x0 = coords[i].x;
				y0 = coords[i].y;
				out.fillCircle(x0, y0, BASE_RADIUS, 1l, getBaseInnerColor(i,
						conf));
				out.setColor(getBaseOuterColor(i, conf));
				out.drawCircle(x0, y0, BASE_RADIUS, 1l);
				out.setColor(getBaseNameColor(i, conf));
				out.drawText(x0, y0, ((ModeleBasesComparison) _listeBases
						.get(i)).getBases());
			}
		} else {
			for (int i = 0; i < _listeBases.size(); i++) {
				ModeleBase mb = _listeBases.get(i);
				x0 = coords[i].x;
				y0 = coords[i].y;
				if (conf._fillBase)
				{
				  out.fillCircle(x0, y0, BASE_RADIUS, 1l, getBaseInnerColor(i,
						conf));
				}
				if (conf._drawOutlineBase)
				{
				out.setColor(getBaseOuterColor(i, conf));
				out.drawCircle(x0, y0, BASE_RADIUS, 1l);
				}
				out.setColor(getBaseNameColor(i, conf));
				out.drawText(x0, y0, _listeBases.get(i).getContent());
			}
		}

		// Drawing base numbers
		double numFontSize = (double) (1.5 * BASE_RADIUS);
		out.setFont(PSExport.FONT_HELVETICA_BOLD, numFontSize);

		for (int i = 0; i < _listeBases.size(); i++) {
			int basenum = _listeBases.get(i).getBaseNumber();
			if (basenum == -1) {
				basenum = i + 1;
			}
			if (this.isNumberDrawn(_listeBases.get(i),conf._numPeriod)) {
				out.setColor(_listeBases.get(i).getStyleBase().get_base_number_color());
				x0 = coords[i].x;
				y0 = coords[i].y;
				x1 = centers[i].x;
				y1 = centers[i].y;
				dx = x1 - x0;
				dy = y1 - y0;
				norm = Math.sqrt(dx * dx + dy * dy);
				dx /= norm;
				dy /= norm;
				out.drawLine((x0 - 1.5 * BASE_RADIUS * dx), (y0 - 1.5
						* BASE_RADIUS * dy), (x0 - 2.5 * BASE_RADIUS * dx),
						(y0 - 2.5 * BASE_RADIUS * dy), 1);
				out.drawText((x0 - (conf._distNumbers+1.0) * BASE_RADIUS * dx), (y0 - (conf._distNumbers+1.0) * BASE_RADIUS
						* dy), "" + (basenum));
			}
		}
		renderAnnotations(out, minX, minY,conf);
		
		// Draw color map 
		if (conf._drawColorMap)
		{ drawColorMap(conf, out); }

		
		// Drawing Title
		Rectangle2D.Double currentBBox = out.getBoundingBox(); 
		double titleFontSize = (2.0*conf._titleFont.getSize());
		out.setColor(conf._titleColor);
		out.setFont(PSExport.FONT_HELVETICA,titleFontSize);
		double yTitle = currentBBox.y-titleFontSize/2.0;
		if (!getName().equals(""))
		{
			out.drawText((maxX-minX)/2.0,yTitle, getName());			
		}
		else if (!conf._title.equals(""))
		{
			out.drawText((maxX-minX)/2.0,yTitle, getName());			
		}


		OutputStreamWriter fout;

		try {
			fout = new OutputStreamWriter(new FileOutputStream(path),"UTF-8");
			
			fout.write(out.export());
			fout.close();
		} catch (IOException e) {
			throw new ExceptionWritingForbidden(e.getMessage());
		}
	}
	
	Point2D.Double buildCaptionPosition(ModeleBase mb, double heightEstimate, VARNAConfig conf)
	{
		double radius = 2.0;
		if (isNumberDrawn(mb,conf._numPeriod))
		{ radius += (conf._distNumbers+1.0); }
		Point2D.Double center = mb.getCenter();
		Point2D.Double p = mb.getCoords();
		double realDistance = BASE_RADIUS*radius + heightEstimate;
		return new Point2D.Double(center.getX()
		+ (p.getX() - center.getX())
		* ((p.distance(center) + realDistance) / p
				.distance(center)), center.getY()
		+ (p.getY() - center.getY())
		* ((p.distance(center) + realDistance) / p
				.distance(center)) );
	}

    public double getBPHeightIncrement()
    {
    	return this._bpHeightIncrement;
    }
	
    public void setBPHeightIncrement(double d)
    {
    	_bpHeightIncrement = d;
    }
	
	public static double CHEM_PROB_ARROW_THICKNESS = 2.0;
	
	private void drawChemProbAnnotation(SecStrDrawingProducer out, ChemProbAnnotation cpa, Point2D.Double anchor, double minX, double minY)
	{
		out.setColor(cpa.getColor());
		Point2D.Double v = cpa.getDirVector();
		Point2D.Double vn = cpa.getNormalVector();
		Point2D.Double base = new Point2D.Double((anchor.x+CHEM_PROB_DIST*v.x),(anchor.y+CHEM_PROB_DIST*v.y));
		Point2D.Double edge = new Point2D.Double((base.x+CHEM_PROB_BASE_LENGTH*cpa.getIntensity()*v.x),(base.y+CHEM_PROB_BASE_LENGTH*cpa.getIntensity()*v.y));
		double thickness = CHEM_PROB_ARROW_THICKNESS*cpa.getIntensity();
		switch (cpa.getType())
		{
		  case ARROW_TYPE:
		  {
			  Point2D.Double arrowTip1 = new Point2D.Double((base.x+cpa.getIntensity()*(CHEM_PROB_ARROW_WIDTH*vn.x+CHEM_PROB_ARROW_HEIGHT*v.x)),
					  (base.y+cpa.getIntensity()*(CHEM_PROB_ARROW_WIDTH*vn.y+CHEM_PROB_ARROW_HEIGHT*v.y)));
			  Point2D.Double arrowTip2 = new Point2D.Double((base.x+cpa.getIntensity()*(-CHEM_PROB_ARROW_WIDTH*vn.x+CHEM_PROB_ARROW_HEIGHT*v.x)),
					  (base.y+cpa.getIntensity()*(-CHEM_PROB_ARROW_WIDTH*vn.y+CHEM_PROB_ARROW_HEIGHT*v.y)));
			  out.drawLine(base.x-minX,minY-base.y,edge.x-minX,minY-edge.y,thickness);
			  out.drawLine(base.x-minX,minY-base.y,arrowTip1.x-minX,minY-arrowTip1.y,thickness);
			  out.drawLine(base.x-minX,minY-base.y,arrowTip2.x-minX,minY-arrowTip2.y,thickness);
		  }
		  break;
		  case PIN_TYPE:
		  {
			  Point2D.Double side1 = new Point2D.Double((edge.x-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.x)),
					  (edge.y-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.y)));
			  Point2D.Double side2 = new Point2D.Double((edge.x-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.x)),
					  (edge.y-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.y)));
			  Point2D.Double side3 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.y)));
			  Point2D.Double side4 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.y)));
				GeneralPath p2 = new GeneralPath();
				p2.moveTo((float)(side1.x-minX),(float)(minY-side1.y));
				p2.lineTo((float)(side2.x-minX),(float)(minY-side2.y));
				p2.lineTo((float)(side3.x-minX),(float)(minY-side3.y));
				p2.lineTo((float)(side4.x-minX),(float)(minY-side4.y));
				p2.closePath();
				out.fillPolygon(p2, cpa.getColor());
				out.drawLine(base.x-minX,minY-base.y,edge.x-minX,minY-edge.y,thickness);
		  }
		  break;
		  case TRIANGLE_TYPE:
		  {
			  Point2D.Double arrowTip1 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_TRIANGLE_WIDTH*vn.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_TRIANGLE_WIDTH*vn.y)));
			  Point2D.Double arrowTip2 = new Point2D.Double((edge.x+cpa.getIntensity()*(-CHEM_PROB_TRIANGLE_WIDTH*vn.x)),
					  (edge.y+cpa.getIntensity()*(-CHEM_PROB_TRIANGLE_WIDTH*vn.y)));
				GeneralPath p2 = new GeneralPath();
				p2.moveTo((float)(base.x-minX),(float)(minY-base.y));
				p2.lineTo((float)(arrowTip1.x-minX),(float)(minY-arrowTip1.y));
				p2.lineTo((float)(arrowTip2.x-minX),(float)(minY-arrowTip2.y));
				p2.closePath();
				out.fillPolygon(p2, cpa.getColor());
		  }
		  break;
		  case DOT_TYPE:
		  {
			  Double radius = CHEM_PROB_DOT_RADIUS*cpa.getIntensity();
			  Point2D.Double center = new Point2D.Double((base.x+radius*v.x)-minX,minY-(base.y+radius*v.y));
			  out.fillCircle(center.x, center.y, radius, thickness, cpa.getColor());
		  }
		  break;
		}
	}

	
	private void renderAnnotations(SecStrDrawingProducer out, double minX, double minY, VARNAConfig conf) 
	{
		for (TextAnnotation textAnnotation : getAnnotations()) {
			out.setColor(textAnnotation.getColor());
			out.setFont(PSExport.FONT_HELVETICA_BOLD,2.0*textAnnotation.getFont().getSize());
			Point2D.Double position = textAnnotation.getCenterPosition();
			if (textAnnotation.getType()==TextAnnotation.BASE)
			{
				ModeleBase mb = (ModeleBase) textAnnotation.getAncrage();
				double fontHeight = Math.ceil(textAnnotation.getFont().getSize());
				position = buildCaptionPosition(mb,fontHeight,conf);

			}
			out.drawText(position.x-minX,-(position.y-minY), textAnnotation.getTexte());
		}
		for (ChemProbAnnotation cpa : getChemProbAnnotations()) {
			Point2D.Double anchor = cpa.getAnchorPosition();
			drawChemProbAnnotation(out,cpa,anchor,minX,minY);
		}
	}


	public boolean isNumberDrawn(ModeleBase mb, int numPeriod)
	{
		if (numPeriod <= 0 )
			return false;
		return  ((mb.getIndex() == 0) || ((mb.getBaseNumber()) % numPeriod == 0)
		|| (mb.getIndex() == get_listeBases().size() - 1));
	}
	
	
	public void saveRNAEPS(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		PSExport out = new PSExport();
		saveRNA(path, conf, 0.4, out);
	}

	public void saveRNAXFIG(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		XFIGExport out = new XFIGExport();
		saveRNA(path, conf, 20, out);
	}

	public void saveRNASVG(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		SVGExport out = new SVGExport();
		saveRNA(path, conf, 0.5, out);
	}

	public Rectangle2D.Double getBBox() {
		Rectangle2D.Double result = new Rectangle2D.Double(10, 10, 10, 10);
		double minx, maxx, miny, maxy;
		minx = Double.MAX_VALUE;
		miny = Double.MAX_VALUE;
		maxx = -Double.MAX_VALUE;
		maxy = -Double.MAX_VALUE;
		for (int i = 0; i < _listeBases.size(); i++) {
			minx = Math.min(_listeBases.get(i).getCoords().getX()
					- BASE_RADIUS, minx);
			miny = Math.min(_listeBases.get(i).getCoords().getY()
					- BASE_RADIUS, miny);
			maxx = Math.max(_listeBases.get(i).getCoords().getX()
					+ BASE_RADIUS, maxx);
			maxy = Math.max(_listeBases.get(i).getCoords().getY()
					+ BASE_RADIUS, maxy);
		}
		result.x = minx;
		result.y = miny;
		result.width = Math.max(maxx - minx,1);
		result.height = Math.max(maxy - miny,1);
		if (_drawMode==RNA.DRAW_MODE_LINEAR)
		{
			double realHeight = _bpHeightIncrement*result.width/2.0;
			result.height += realHeight;
			result.y -= realHeight; 
		}
		return result;
	}

	public void setCoord(int index, Point2D.Double p) {
		setCoord(index, p.x, p.y);
	}

	public void setCoord(int index, double x, double y) {
		if (index < _listeBases.size()) {
			_listeBases.get(index).setCoords(new Point2D.Double(x, y));
		}
	}

	public Point2D.Double getCoords(int i) {
		if (i < _listeBases.size() && i >= 0) {
			return _listeBases.get(i).getCoords();
		}
		return new Point2D.Double();
	}

	public String getBaseContent(int i) {
		if ((i>=0)&&(i < _listeBases.size())) 
		{ return _listeBases.get(i).getContent(); }
		return "";
	}

	public int getBaseNumber(int i) {
		if ((i>=0)&&(i < _listeBases.size())) 
		{ return _listeBases.get(i).getBaseNumber(); }
		return -1;
	}
	
	public Point2D.Double getCenter(int i) {
		if (i < _listeBases.size()) {
			return _listeBases.get(i).getCenter();
		}

		return new Point2D.Double();
	}

	public void setCenter(int i, double x, double y) {
		setCenter(i, new Point2D.Double(x,y));
	}

	public void setCenter(int i, Point2D.Double p) {
		if (i < _listeBases.size()) {
			_listeBases.get(i).setCenter(p);
		}
	}

	public void drawRNACircle() {
		_drawn = true;
		_drawMode = DRAW_MODE_CIRCULAR;
		int radius = (int) ((3 * (_listeBases.size() + 1) * BASE_RADIUS) / (2 * Math.PI));
		double angle;
		for (int i = 0; i < _listeBases.size(); i++) {
			angle = -((((double) -(i + 1)) * 2.0 * Math.PI)
					/ ((double) (_listeBases.size() + 1)) - Math.PI / 2.0);
			_listeBases.get(i).setCoords(
					new Point2D.Double(
							(radius * Math.cos(angle) * _spaceBetweenBases),
							(radius * Math.sin(angle) * _spaceBetweenBases)));
			_listeBases.get(i).setCenter(new Point2D.Double(0, 0));
		}
	}

	public void drawRNAVARNAView() {
		_drawn = true;
		_drawMode = DRAW_MODE_VARNA_VIEW;
		VARNASecDraw vs = new VARNASecDraw();
		vs.drawRNA(1, this);
	}
	
	
	public void drawRNALine() {
		_drawn = true;
		_drawMode = DRAW_MODE_LINEAR;
		for (int i = 0; i < get_listeBases().size(); i++) {
			get_listeBases().get(i).setCoords(
					new Point2D.Double(i * _spaceBetweenBases * 20, 0));
			get_listeBases().get(i).setCenter(
					new Point2D.Double(i * _spaceBetweenBases * 20, -10));
		}
	}
	
	
	/**
	 * IN: Argument helixEndPoint is an IN argument (will be read),
	 * and must contain an helix edge endpoint.
	 * 
	 * The other arguments are OUT arguments
	 * (must be existing objects, content will be overwritten).
	 * 
	 * OUT: The i argument will contain a vector colinear to the vector
	 * from the helix startPosition to endPosition or the opposite
	 * depending on there the endpoint is (the endpoint will be on the
	 * destination side of the vector). ||i|| = 1
	 * 
	 * OUT: The j vector will contain an vector that is colinear
	 * to the last/first base pair connection on the side of this endpoint.
	 * The vector will be oriented to the side of the given endpoint.
	 * ||j|| = 1
	 */
	private void computeHelixEndPointDirections(
			EdgeEndPoint helixEndPoint, // IN
			Point2D.Double i, // OUT
			Point2D.Double j // OUT
			) {
		RNATemplateHelix helix = (RNATemplateHelix) helixEndPoint.getElement();
		Point2D.Double startpos = helix.getStartPosition();
		Point2D.Double endpos = helix.getEndPosition();
		Point2D.Double helixVector = new Point2D.Double();
		switch (helixEndPoint.getPosition()) {
		case IN1:
		case OUT2:
			helixVector.x = startpos.x - endpos.x;
			helixVector.y = startpos.y - endpos.y;
			break;
		case IN2:
		case OUT1:
			helixVector.x = endpos.x - startpos.x;
			helixVector.y = endpos.y - startpos.y;
			break;
		}
		double helixVectorLength = Math.hypot(helixVector.x, helixVector.y);
		// i is the vector which is colinear to helixVector and such that ||i|| = 1
		i.x = helixVector.x / helixVectorLength;
		i.y = helixVector.y / helixVectorLength;
		// Find j such that it is orthogonal to i, ||j|| = 1
		// and j goes to the side where the sequence will be connected
		switch (helixEndPoint.getPosition()) {
		case IN1:
		case IN2:
			// rotation of +pi/2
			j.x = - i.y;
			j.y =   i.x;
			break;
		case OUT1:
		case OUT2:
			// rotation of -pi/2
			j.x =   i.y;
			j.y = - i.x;
			break;
		}
		if (helix.isFlipped()) {
			j.x = - j.x;
			j.y = - j.y;
		}

	}
	
	/**
	 * A cubic Bezier curve can be defined by 4 points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * For each of the curve end points, there is the last/first point of the
	 * curve and a point which gives the direction and length of the tangent
	 * vector on that side. This two points are respectively curveEndPoint
	 * and curveVectorOtherPoint.
	 * IN:  Argument helixVector is the vector formed by the helix,
	 *      in the right direction for our sequence.
	 * IN:  Argument curveEndPoint is the position of the endpoint on the helix.
	 * OUT: Argument curveVectorOtherPoint must be allocated
	 *      and the values will be modified.
	 */
	private void computeBezierTangentVectorTarget(
			EdgeEndPoint endPoint,
			Point2D.Double curveEndPoint,
			Point2D.Double curveVectorOtherPoint)
			throws RNATemplateDrawingAlgorithmException {
		
		boolean sequenceEndPointIsIn;
		RNATemplateUnpairedSequence sequence;
		
		if (endPoint.getElement() instanceof RNATemplateHelix) {
			sequence = (RNATemplateUnpairedSequence) endPoint.getOtherElement();
			EdgeEndPointPosition endPointPositionOnHelix = endPoint.getPosition();
			switch (endPointPositionOnHelix) {
			case IN1:
			case IN2:
				sequenceEndPointIsIn = false;
				break;
			default:
				sequenceEndPointIsIn = true;
			}
			
			EdgeEndPoint endPointOnHelix =
				sequenceEndPointIsIn ?
						sequence.getIn().getOtherEndPoint() :
						sequence.getOut().getOtherEndPoint();
			if (endPointOnHelix == null) {
				throw (new RNATemplateDrawingAlgorithmException("Sequence is not connected to an helix."));
			}
		} else {
			// The endpoint is on an unpaired sequence.
			sequence = (RNATemplateUnpairedSequence) endPoint.getElement();
			if (endPoint == sequence.getIn()) {
				// endpoint is 5'
				sequenceEndPointIsIn = true;
			} else {
				sequenceEndPointIsIn = false;
			}
		}

		double l =
			sequenceEndPointIsIn ?
				sequence.getInTangentVectorLength() :
				sequence.getOutTangentVectorLength();
		
		// Compute the absolute angle our line makes to the helix
		double theta =
			sequenceEndPointIsIn ?
				sequence.getInTangentVectorAngle() :
				sequence.getOutTangentVectorAngle();
		
		// Compute v, the tangent vector of the Bezier curve
		Point2D.Double v = new Point2D.Double();
		v.x = l * Math.cos(theta);
		v.y = l * Math.sin(theta);
		curveVectorOtherPoint.x = curveEndPoint.x + v.x;
		curveVectorOtherPoint.y = curveEndPoint.y + v.y;
	}
	

	/**
	 * Compute the angle made by a vector.
	 */
	private static double angleFromVector(Point2D.Double v) {
		return angleFromVector(v.x, v.y);
	}
	
	private static double angleFromVector(double x, double y) {
		double l = Math.hypot(x, y);
		if (y > 0) {
			return Math.acos(x / l);
		} else if (y < 0) {
			return - Math.acos(x / l);
		} else {
			return x > 0 ? 0 : Math.PI;
		}
	}
	
	/**
	 * Compute (actual helix length / helix length in template).
	 */
	private double computeLengthIncreaseFactor(
			int[] basesInHelixArray,  // IN
			RNATemplateHelix helix    // IN
			) {
		double templateLength = computeHelixTemplateLength(helix);
		double realLength = computeHelixRealLength(basesInHelixArray);
		return realLength / templateLength;
	}
	
	/**
	 * Compute (actual helix vector - helix vector in template).
	 */
	private Point2D.Double computeLengthIncreaseDelta(
			int[] basesInHelixArray,  // IN
			RNATemplateHelix helix    // IN
			) {
		double templateLength = computeHelixTemplateLength(helix);
		double realLength = computeHelixRealLength(basesInHelixArray);
		Point2D.Double i = new Point2D.Double();
		computeTemplateHelixVectors(helix, null, i, null);
		return new Point2D.Double(i.x*(realLength-templateLength), i.y*(realLength-templateLength));
	}
	
	/**
	 * Compute helix interesting vectors from template helix.
	 * @param helix The template helix you want to compute the vectors from.
	 * @param o This point coordinates will be set the origin of the helix (or not if null),
	 *          ie. the point in the middle of the base pair with the two most extreme bases.
	 * @param i Will be set to the normalized helix vector. (nothing done if null)
	 * @param j Will be set to the normalized helix base pair vector (5' -> 3'). (nothing done if null)
	 */
	private void computeTemplateHelixVectors(
			RNATemplateHelix helix,  // IN
			Point2D.Double o,        // OUT
			Point2D.Double i,        // OUT
			Point2D.Double j         // OUT
			) {
		Point2D.Double startpos, endpos;
		if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
			startpos = helix.getStartPosition();
			endpos = helix.getEndPosition();
		} else {
			endpos = helix.getStartPosition();
			startpos = helix.getEndPosition();
		}
		if (o != null) {
			o.x = startpos.x;
			o.y = startpos.y;
		}
		if (i != null || j != null) {
			// (i_x,i_y) is the vector between two consecutive bases of the same side of an helix
			if (i == null)
				i = new Point2D.Double();
			i.x = (endpos.x - startpos.x);
			i.y = (endpos.y - startpos.y);
			double i_original_norm = Math.hypot(i.x, i.y);
			// change its norm to 1
			i.x = i.x / i_original_norm;
			i.y = i.y / i_original_norm;
			if (j != null) {
				j.x = - i.y;
				j.y =   i.x;
				if (helix.isFlipped()) {
					j.x = - j.x;
					j.y = - j.y;
				}
				double j_original_norm = Math.hypot(j.x, j.y);
				// change (j_x,j_y) so that its norm is 1
				j.x = j.x / j_original_norm;
				j.y = j.y / j_original_norm;
			}
		}
	}
	
	
	private static class ComputeArcCenter {
		
		/**
		 * Given an arc length (l) and segment length (delta) of the arc,
		 * find where to put the center, returned as a position of the perpendicular
		 * bisector of the segment. The positive side is the one where the arc is drawn.
		 * It works using Newton's method.
		 */
		public static double computeArcCenter(double delta, double l) {
			double x_n = 0;
			double x_n_plus_1, f_x_n, f_x_n_plus_1;
			int steps = 0;
			while (true) {
				f_x_n = f(x_n,delta);
				x_n_plus_1 = x_n - (f_x_n - l)/fprime(x_n,delta);
				f_x_n_plus_1 = f(x_n_plus_1,delta);
				steps++;
				// We want a precision of 0.1 on arc length
				if (x_n_plus_1 == Double.NEGATIVE_INFINITY || Math.abs(f_x_n_plus_1 - f_x_n) < 0.1) {
					//System.out.println("computeArcCenter: steps = " + steps + "    result = " + x_n_plus_1);
					return x_n_plus_1;
				}
				x_n = x_n_plus_1;
				f_x_n = f_x_n_plus_1;
			}
		}
		
		public static double f(double c, double delta) {
			if (c < 0) {
				return 2*Math.atan(delta/(-2*c)) * Math.sqrt(delta*delta/4 + c*c);
			} else if (c != 0) { // c > 0
				return (2*Math.PI - 2*Math.atan(delta/(2*c))) * Math.sqrt(delta*delta/4 + c*c);
			} else { // c == 0
				return Math.PI * Math.sqrt(delta*delta/4 + c*c);
			}
		}
		
		/**
		 * d/dc f(c,delta)
		 */
		public static double fprime(double c, double delta) {
			if (c < 0) {
				return delta/(c*c + delta/4)*Math.sqrt(delta*delta/4 + c*c) + 2*Math.atan(delta/(-2*c))*c/Math.sqrt(delta*delta/4 + c*c);
			} else if (c != 0) { // c > 0
				return delta/(c*c + delta/4)*Math.sqrt(delta*delta/4 + c*c) + (2*Math.PI - 2*Math.atan(delta/(-2*c)))*c/Math.sqrt(delta*delta/4 + c*c);
			} else { // c == 0
				return 2;
			}
		}
	}
	
	
	/**
	 * Estimate bulge arc length.
	 */
	private double estimateBulgeArcLength(int firstBase, int lastBase) {
		if (firstBase + 1 == lastBase)
			return LOOP_DISTANCE; // there is actually no bulge
		double len = 0.0;
		int k = firstBase;
		while (k < lastBase) {
			int l = _listeBases.get(k).getElementStructure();
			if (k < l && l < lastBase) {
				len += BASE_PAIR_DISTANCE;
				k = l;
			} else {
				len += LOOP_DISTANCE;
				k++;
			}
		}
		return len;
	}
	
	
	/**
	 * Estimate bulge width, the given first and last bases must be those in the helix.
	 */
	private double estimateBulgeWidth(int firstBase, int lastBase) {
		double len = estimateBulgeArcLength(firstBase, lastBase);
		return 2 * (len / Math.PI);
	}
	
	
	/**
	 * Estimate bulge width, the given first and last bases must be those in the helix.
	 */
	private double estimateBulgeWidthOLD(int firstBase, int lastBase) {
		if (firstBase + 1 == lastBase)
			return LOOP_DISTANCE; // there is actually no bulge
		double len = estimateBulgeArcLength(firstBase, lastBase);
		// We expect the bases to be drawn on a circle having its center in the middle of the helix.
		len += BASE_PAIR_DISTANCE; // 2 * half base pair part of the helix
		return 2 * (len / Math.PI);
	}
	
	
	/**
	 * Get helix length in template.
	 */
	private double computeHelixTemplateLength(RNATemplateHelix helix) {
		return Math.hypot(helix.getStartPosition().x - helix.getEndPosition().x,
				helix.getStartPosition().y - helix.getEndPosition().y);
	}
	
	
	/**
	 * Compute helix actual length (as drawHelixLikeTemplateHelix() would draw it).
	 */
	private double computeHelixRealLength(int[] basesInHelixArray) {
		return drawHelixLikeTemplateHelix(basesInHelixArray, null, null, null, 0, null);
	}
	
	
	/**
	 * Draw the given helix (given as a *SORTED* array of indexes)
	 * like defined in the given template helix.
	 * OUT: The bases positions are not changed in fact,
	 *      instead the coords and centers arrays are modified.
	 * IN:  The helix origin position is multiplied by scaleHelixOrigin
	 *      and translateVectors.get(helix) is added.
	 * RETURN VALUE:
	 *      The length of the drawn helix.
	 * 
	 */
	private double drawHelixLikeTemplateHelix(
			int[] basesInHelixArray,  // IN
			RNATemplateHelix helix,   // IN  (optional, ie. may be null)
			Point2D.Double[] coords,  // OUT (optional, ie. may be null)
			Point2D.Double[] centers, // OUT (optional, ie. may be null)
			double scaleHelixOrigin,  // IN
			Map<RNATemplateHelix, Point2D.Double> translateVectors // IN (optional, ie. may be null)
			) {
		int n = basesInHelixArray.length / 2;
		if (n == 0)
			return 0;
		 // Default values when not template helix is provided:
		Point2D.Double o = new Point2D.Double(0, 0);
		Point2D.Double i = new Point2D.Double(1, 0);
		Point2D.Double j = new Point2D.Double(0, 1);
		boolean flipped = false;
		if (helix != null) {
			computeTemplateHelixVectors(helix, o, i, j);
			flipped = helix.isFlipped();
		}
		Point2D.Double li = new Point2D.Double(i.x*LOOP_DISTANCE, i.y*LOOP_DISTANCE);
		// We want o to be the point where the first base (5' end) is
		o.x = (o.x - j.x * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
		o.y = (o.y - j.y * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
		if (translateVectors != null && translateVectors.containsKey(helix)) {
			Point2D.Double v = translateVectors.get(helix);
			o.x = o.x + v.x;
			o.y = o.y + v.y;
		}
		
		// We need this array so that we can store positions even if coords == null
		Point2D.Double[] helixBasesPositions = new Point2D.Double[basesInHelixArray.length];
		for (int k=0; k<helixBasesPositions.length; k++) {
			helixBasesPositions[k] = new Point2D.Double();	
		}
		Point2D.Double accDelta = new Point2D.Double(0, 0);
		for (int k=0; k<n; k++) {
			int kp = 2*n-k-1;
			Point2D.Double p1 = helixBasesPositions[k]; // we assign the point *reference*
			Point2D.Double p2 = helixBasesPositions[kp];
			// Do we have a bulge between previous base pair and this one?
			boolean bulge = k >= 1 && (basesInHelixArray[k] != basesInHelixArray[k-1] + 1
					                   || basesInHelixArray[kp+1] != basesInHelixArray[kp] + 1);
			if (k >= 1) {
				if (basesInHelixArray[k] < basesInHelixArray[k-1]
				    || basesInHelixArray[kp+1] < basesInHelixArray[kp]) {
					throw new Error("Internal bug: basesInHelixArray must be sorted");
				}
				if (bulge) {
					// Estimate a good distance (delta) between the previous base pair and this one
					double delta1 = estimateBulgeWidth(basesInHelixArray[k-1], basesInHelixArray[k]);
					double delta2 = estimateBulgeWidth(basesInHelixArray[kp], basesInHelixArray[kp+1]);
					// The biggest bulge defines the width
					double delta = Math.max(delta1, delta2);

					if (coords != null) {
						// Now, where do we put the bases that are part of the bulge?
						for (int side=0; side<2; side++) {
							Point2D.Double pstart = new Point2D.Double();
							Point2D.Double pend = new Point2D.Double();
							Point2D.Double bisectVect = new Point2D.Double();
							Point2D.Double is = new Point2D.Double();
							int firstBase, lastBase;
							double alphasign = flipped ? -1 : 1;
							if (side == 0) {
								firstBase = basesInHelixArray[k-1];
								lastBase  = basesInHelixArray[k];
								pstart.setLocation(o.x + accDelta.x,
								                   o.y + accDelta.y);
								pend.setLocation(o.x + accDelta.x + i.x*delta,
								                 o.y + accDelta.y + i.y*delta);
								bisectVect.setLocation(-j.x, -j.y);
								is.setLocation(i);
							} else {
								firstBase = basesInHelixArray[kp];
								lastBase  = basesInHelixArray[kp+1];
								pstart.setLocation(o.x + accDelta.x + i.x*delta + j.x*BASE_PAIR_DISTANCE,
						                           o.y + accDelta.y + i.y*delta + j.y*BASE_PAIR_DISTANCE);
								pend.setLocation(o.x + accDelta.x + j.x*BASE_PAIR_DISTANCE,
								                 o.y + accDelta.y + j.y*BASE_PAIR_DISTANCE);

								bisectVect.setLocation(j);
								is.setLocation(-i.x, -i.y);
							}
							double arclen = estimateBulgeArcLength(firstBase, lastBase);
							double centerOnBisect = ComputeArcCenter.computeArcCenter(delta, arclen);
							
							// Should we draw the base on an arc or simply use a line?
							if (centerOnBisect > -1000) {
								Point2D.Double center = new Point2D.Double(pstart.x + is.x*delta/2 + bisectVect.x*centerOnBisect,
								                                           pstart.y + is.y*delta/2 + bisectVect.y*centerOnBisect);
								int b = firstBase;
								double len = 0;
								double r = Math.hypot(pstart.x - center.x, pstart.y - center.y);
								double alpha0 = angleFromVector(pstart.x - center.x, pstart.y - center.y);
								while (b < lastBase) {
									int l = _listeBases.get(b).getElementStructure();
									if (b < l && l < lastBase) {
										len += BASE_PAIR_DISTANCE;
										b = l;
									} else {
										len += LOOP_DISTANCE;
										b++;
									}
									if (b < lastBase) {
										coords[b].x = center.x + r*Math.cos(alpha0 + alphasign*len/r);
										coords[b].y = center.y + r*Math.sin(alpha0 + alphasign*len/r);
									}
								}
							} else {
								// Draw on a line
								double nBP = 0;
								double nLD = 0;
								{
									int b = firstBase;
									while (b < lastBase) {
										int l = _listeBases.get(b).getElementStructure();
										if (b < l && l < lastBase) {
											nBP++;
											b = l;
										} else {
											nLD++;
											b++;
										}
									}
								}
								// Distance between paired bases cannot be changed
								// (imposed by helix width) but distance between other
								// bases can be adjusted.
								double LD = Math.max((delta - nBP*BASE_PAIR_DISTANCE) / nLD, 0);
								//System.out.println("nBP=" + nBP + " nLD=" + nLD);
								double len = 0;
								{
									int b = firstBase;
									while (b < lastBase) {
										int l = _listeBases.get(b).getElementStructure();
										if (b < l && l < lastBase) {
											len += BASE_PAIR_DISTANCE;
											b = l;
										} else {
											len += LD;
											b++;
										}
										if (b < lastBase) {
											coords[b].x = pstart.x + is.x*len;
											coords[b].y = pstart.y + is.y*len;
										}
									}
								}
								//System.out.println("len=" + len + " delta=" + delta + " d(pstart,pend)=" + Math.hypot(pend.x-pstart.x, pend.y-pstart.y));
							}
							
							// Does the bulge contain an helix?
							// If so, use drawLoop() to draw it.
							{
								int b = firstBase;
								while (b < lastBase) {
									int l = _listeBases.get(b).getElementStructure();
									if (b < l && l < lastBase) {
										// Helix present in bulge
										Point2D.Double b1pos = coords[b];
										Point2D.Double b2pos = coords[l];
										double beta = angleFromVector(b2pos.x - b1pos.x, b2pos.y - b1pos.y) - Math.PI / 2 + (flipped ? Math.PI : 0);
										Point2D.Double loopCenter = new Point2D.Double((b1pos.x + b2pos.x)/2, (b1pos.y + b2pos.y)/2);
										drawLoop(b,
												 l,
												 loopCenter.x,
												 loopCenter.y,
												 beta,
												 coords,
												 centers);
										// If the helix is flipped, we need to compute the symmetric
										// of the whole loop.
										if (helix.isFlipped()) {
											Point2D.Double v = new Point2D.Double(Math.cos(beta), Math.sin(beta));
											Point2D.Double[] points1 = new Point2D.Double[l-b+1];
											Point2D.Double[] points2 = new Point2D.Double[l-b+1];
											for (int c=b; c<=l; c++) {
												// This is an assignment by reference.
												points1[c-b] = coords[c];
												points2[c-b] = centers[c];
											}
											symmetric(loopCenter, v, points1);
											symmetric(loopCenter, v, points2);
										}
										// Continue
										b = l;
									} else {
										b++;
									}
								}
							}
						}
					}
					
					accDelta.x += i.x*delta;
					accDelta.y += i.y*delta;
					p1.x = o.x + accDelta.x;
					p1.y = o.y + accDelta.y;
					p2.x = p1.x + j.x*BASE_PAIR_DISTANCE;
					p2.y = p1.y + j.y*BASE_PAIR_DISTANCE;
					
				} else {
					accDelta.x += li.x;
					accDelta.y += li.y;
					p1.x = o.x + accDelta.x;
					p1.y = o.y + accDelta.y;
					p2.x = p1.x + j.x*BASE_PAIR_DISTANCE;
					p2.y = p1.y + j.y*BASE_PAIR_DISTANCE;
				}
			} else {
				// First base pair
				p1.x = o.x;
				p1.y = o.y;
				p2.x = p1.x + j.x*BASE_PAIR_DISTANCE;
				p2.y = p1.y + j.y*BASE_PAIR_DISTANCE;
			}
		}
		
		Point2D.Double p1 = helixBasesPositions[0];
		Point2D.Double p2 = helixBasesPositions[n-1];
		
		if (coords != null) {
			for (int k=0; k<helixBasesPositions.length; k++) {
				coords[basesInHelixArray[k]] = helixBasesPositions[k];
			}
		}
		
		return Math.hypot(p2.x-p1.x, p2.y-p1.y);
	}
	

	
	private void drawHelixLikeTemplateHelixOLD(
			int[] basesInHelixArray,  // IN
			RNATemplateHelix helix,   // IN
			Point2D.Double[] coords,  // OUT
			Point2D.Double[] centers, // OUT
			double scaleHelixOrigin,  // IN
			Map<RNATemplateHelix, Point2D.Double> translateVectors // IN
			) {
		int n = basesInHelixArray.length / 2;
		Point2D.Double o = new Point2D.Double();
		Point2D.Double i = new Point2D.Double();
		Point2D.Double j = new Point2D.Double();
		computeTemplateHelixVectors(helix, o, i, j);
		i.x = i.x * LOOP_DISTANCE;
		i.y = i.y * LOOP_DISTANCE;
		o.x = (o.x - j.x * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
		o.y = (o.y - j.y * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
		if (translateVectors != null && translateVectors.containsKey(helix)) {
			Point2D.Double v = translateVectors.get(helix);
			o.x = o.x + v.x;
			o.y = o.y + v.y;
		}
		
		for (int k=0; k<n; k++) {
			int b1 = basesInHelixArray[k];
			int b2 = basesInHelixArray[2*n-k-1];
			Point2D.Double p1 = coords[b1]; // we assign the point *reference*
			p1.x = o.x + k*i.x;
			p1.y = o.y + k*i.y;
			Point2D.Double p2 = coords[b2];
			p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
			p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;
		}
	
		for (int k=0; k<2*n-1; k++) {
			if (k == n-1) continue;
			int b1 = basesInHelixArray[k];
			int b2 = basesInHelixArray[k+1];
			if (b1 + 1 != b2) {
				// There is a loop between these 2 bases
				Point2D.Double b1pos = coords[b1];
				Point2D.Double b2pos = coords[b2];
				// ||v|| = 1
				Point2D.Double v = new Point2D.Double();
				if (k >= n) {
					v.x = j.x;
					v.y = j.y;
				} else {
					v.x = - j.x;
					v.y = - j.y;
				}
				Point2D.Double loopCenter = new Point2D.Double();
				loopCenter.x = (b1pos.x + b2pos.x)/2 + v.x * LOOP_DISTANCE;
				loopCenter.y = (b1pos.y + b2pos.y)/2 + v.y * LOOP_DISTANCE;
				drawLoop(b1+1,
						 b2-1,
						 loopCenter.x,
						 loopCenter.y,
						 angleFromVector(v),
						 coords,
						 centers);
				// If the helix is flipped, we need to compute the symmetric
				// of the whole loop.
				if (helix.isFlipped()) {
					int from = b1+1;
					int to = b2-1;
					Point2D.Double[] points1 = new Point2D.Double[to-from+1];
					Point2D.Double[] points2 = new Point2D.Double[to-from+1];
					for (int b=from; b<=to; b++) {
						// This is an assignment by reference.
						points1[b-from] = coords[b];
						points2[b-from] = centers[b];
					}
					symmetric(loopCenter, v, points1);
					symmetric(loopCenter, v, points2);
				}
			}
		}
	}
	
	/**
	 * A Bezier curve can be defined by four points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * Here we give this four points and an array of bases indexes
	 * (which must be indexes in this RNA sequence) which will be moved
	 * to be on the Bezier curve.
	 * The bases positions are not changed in fact, instead the coords and
	 * centers arrays are modified.
	 */
	private void drawOnBezierCurve(int[] basesInSequence,
			Point2D.Double P0,
			Point2D.Double P1,
			Point2D.Double P2,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		// Draw the bases of the sequence along a Bezier curve
		int n = basesInSequence.length;
		// We choose to approximate the Bezier curve by 10*n straight lines.
		CubicBezierCurve bezier = new CubicBezierCurve(P0, P1, P2, P3, 10*n);
		double curveLength = bezier.getApproxCurveLength();
		double delta_t = curveLength / (n+1);
		double[] t = new double[n];
		for (int k=0; k<n; k++) {
			t[k] = (k+1) * delta_t;
		}
		Point2D.Double[] sequenceBasesCoords = bezier.uniformParam(t);
		for (int k=0; k<n; k++) {
			coords[basesInSequence[k]] = sequenceBasesCoords[k];
		}
	}
	
	/**
	 * Like drawOnBezierCurve(), but on a straight line.
	 */
	private void drawOnStraightLine(int[] basesInSequence,
			Point2D.Double P0,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		// Draw the bases of the sequence along a Bezier curve
		int n = basesInSequence.length;
		Point2D.Double v = new Point2D.Double(P3.x - P0.x, P3.y - P0.y);
		for (int k=0; k<n; k++) {
			coords[basesInSequence[k]].x = P0.x + (k+1) / (double) (n+1) * v.x;
			coords[basesInSequence[k]].y = P0.y + (k+1) / (double) (n+1) * v.y;
		}
	}
	
	/**
	 * This functions draws the RNA sequence between (including)
	 * firstBase and lastBase along a curve.
	 * The sequence may contain helices.
	 * 
	 * Bezier curve:
	 * A Bezier curve can be defined by four points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * 
	 * Straight line:
	 * If P1 and P2 are null, the bases are drawn on a straight line.
	 * 
	 * OUT: The bases positions are not changed in fact,
	 * instead the coords and centers arrays are modified.
	 * 
	 */
	private void drawAlongCurve(
			int firstBase,
			int lastBase,
			Point2D.Double P0,
			Point2D.Double P1,
			Point2D.Double P2,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		
		// First we find the bases which are directly on the Bezier curve
		ArrayList<Integer> alongBezierCurve = new ArrayList<Integer>();
		for (int depth=0, i=firstBase; i<=lastBase; i++) {
			int k = _listeBases.get(i).getElementStructure();
			if (k < 0 || k > lastBase || k < firstBase) {
				if (depth == 0) {
					alongBezierCurve.add(i);
				}
			} else {
				if (i < k) {
					if (depth == 0) {
						alongBezierCurve.add(i);
						alongBezierCurve.add(k);
					}
					depth++;
				} else {
					depth--;
				}
			}
		}
		// number of bases along the Bezier curve
		int n = alongBezierCurve.size();
		int[] alongBezierCurveArray = RNATemplateAlign.intArrayFromList(alongBezierCurve);
		if (n > 0) {
			if (P1 != null && P2 != null) {
				drawOnBezierCurve(alongBezierCurveArray, P0, P1, P2, P3, coords, centers);
			} else {
				drawOnStraightLine(alongBezierCurveArray, P0, P3, coords, centers);
			}
		}
		// Now use the radiate algorithm to draw the helixes
		for (int k=0; k<n-1; k++) {
			int b1 = alongBezierCurveArray[k];
			int b2 = alongBezierCurveArray[k+1];
			if (_listeBases.get(b1).getElementStructure() == b2) {
				Point2D.Double b1pos = coords[b1];
				Point2D.Double b2pos = coords[b2];
				double alpha = angleFromVector(b2pos.x - b1pos.x, b2pos.y - b1pos.y);
				drawLoop(b1,
						 b2,
						 (b1pos.x + b2pos.x)/2,
						 (b1pos.y + b2pos.y)/2,
						 alpha - Math.PI / 2,
						 coords,
						 centers);
			}
		}	
	}
	
	
	/**
	 * Compute the symmetric of all the points in the points array
	 * relative to the line that goes through p and has director vector v.
	 * The array is modified in-place.
	 */
	private static void symmetric(
			Point2D.Double p,
			Point2D.Double v,
			Point2D.Double[] points) {
		// ||v||^2
		double lv = v.x*v.x + v.y*v.y;
		for (int i=0; i<points.length; i++) {
			// A is the coordinates of points[i] after moving the origin at p
			Point2D.Double A = new Point2D.Double(points[i].x - p.x, points[i].y - p.y);
			// Symmetric of A
			Point2D.Double B = new Point2D.Double(
					-(A.x*v.y*v.y - 2*A.y*v.x*v.y - A.x*v.x*v.x) / lv,
					 (A.y*v.y*v.y + 2*A.x*v.x*v.y - A.y*v.x*v.x) / lv);
			// Change the origin back
			points[i].x = B.x + p.x;
			points[i].y = B.y + p.y;
		}
	}
	
	private void computeHelixTranslations(
			Tree<RNANodeValueTemplate> tree, // IN
			Map<RNATemplateHelix, Point2D.Double> translateVectors, // OUT (must be initialized)
			RNATemplateMapping mapping,      // IN
			Point2D.Double parentDeltaVector // IN
	) {
		RNANodeValueTemplate nvt = tree.getValue();
		Point2D.Double newDeltaVector = parentDeltaVector;
		if (nvt instanceof RNANodeValueTemplateBasePair) {
			RNATemplateHelix helix = ((RNANodeValueTemplateBasePair) nvt).getHelix();
			if (! translateVectors.containsKey(helix)) {
				translateVectors.put(helix, parentDeltaVector);
				int[] basesInHelixArray;
				if (mapping.getAncestor(helix) != null) {
					basesInHelixArray = RNATemplateAlign.intArrayFromList(mapping.getAncestor(helix));
					Arrays.sort(basesInHelixArray);
				} else {
					basesInHelixArray = new int[0];
				}
				Point2D.Double helixDeltaVector = computeLengthIncreaseDelta(basesInHelixArray, helix);
				newDeltaVector = new Point2D.Double(parentDeltaVector.x+helixDeltaVector.x, parentDeltaVector.y+helixDeltaVector.y);
			} 
		}
		for (Tree<RNANodeValueTemplate> subtree: tree.getChildren()) {
			computeHelixTranslations(subtree, translateVectors, mapping, newDeltaVector);
		}
	}
	
	private Map<RNATemplateHelix, Point2D.Double> computeHelixTranslations(
			Tree<RNANodeValueTemplate> tree, // IN
			RNATemplateMapping mapping       // IN
	) {
		Map<RNATemplateHelix, Point2D.Double> translateVectors = new HashMap<RNATemplateHelix, Point2D.Double>();
		computeHelixTranslations(tree, translateVectors, mapping, new Point2D.Double(0,0));
		return translateVectors;
	}
	
	
	/**
	 * Same as below, with default helixLengthAdjustmentMethod.
	 */
	public RNATemplateMapping drawRNATemplate(RNATemplate template) throws RNATemplateDrawingAlgorithmException {
		return drawRNATemplate(template, DrawRNATemplateMethod.NOADJUST);
	}
	
	
	/**
	 * What to do in case some helices are of a different length
	 * in the template and the actual helix.
	 * Possibles values are:
	 * 0 - No adjustment is done.
	 *     A longer than expected helix might bump into an other helix.
	 * 1 - Scaling factors (actual length / template length) are calculated,
	 *     the maximum scaling factor L is extracted, then all helix positions
	 *     are multiplied by L. 
	 * 2 - Same as 1, but L is computed as the minimum value such that there
	 *     are no backbone intersections.
	 */
	public static class DrawRNATemplateMethod {
		public static int NOADJUST = 0;
		public static int MAXSCALINGFACTOR = 1;
		public static int NOINTERSECT = 2;
		public static int HELIXTRANSLATE = 3;
	}
	
	/**
	 * Draw this RNA like the given template.
	 * The helixLengthAdjustmentMethod argument tells what to do in case
	 * some helices are of a different length in the template and the
	 * actual helix. See class DrawRNATemplateMethod above for possible values.
	 */
	public RNATemplateMapping drawRNATemplate(
			RNATemplate template,
			int helixLengthAdjustmentMethod)
	throws RNATemplateDrawingAlgorithmException {
		_drawn = true;
		_drawMode = DRAW_MODE_TEMPLATE;
		
		// debug
//		try {
//			RNA perfectMatchingRNA = template.toRNA();
//			System.out.println("An RNA that would perfectly match this template would be:");
//			System.out.println(perfectMatchingRNA.getStructDBN());
//		} catch (ExceptionInvalidRNATemplate e) {
//			e.printStackTrace();
//		}
		
		RNATemplateMapping mapping = RNATemplateAlign.mapRNAWithTemplate(this, template);
		System.out.println(mapping.showCompact(this));
		
		// debug
//			RNATemplateAlign.printMapping(mapping, template, getSeq());
//			try {
//				TreeGraphviz.treeToGraphvizPostscript(alignment, "alignment_graphviz.ps");
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
		
		

		Iterator<RNATemplateElement> iter;
		double globalIncreaseFactor = 1;
		Map<RNATemplateHelix, Point2D.Double> translateVectors = null;
		if (helixLengthAdjustmentMethod == DrawRNATemplateMethod.MAXSCALINGFACTOR) {
			// Compute increase factors for helices.
			Map<RNATemplateHelix, Double> lengthIncreaseFactor = new HashMap<RNATemplateHelix, Double>();
			double maxLengthIncreaseFactor = Double.NEGATIVE_INFINITY;
			RNATemplateHelix maxIncreaseHelix = null;
			iter = template.rnaIterator();
			while (iter.hasNext()) {
				RNATemplateElement element = iter.next();
				if (element instanceof RNATemplateHelix
						&& mapping.getAncestor(element) != null
						&& !lengthIncreaseFactor.containsKey(element)) {
					RNATemplateHelix helix = (RNATemplateHelix) element;
					int[] basesInHelixArray = RNATemplateAlign.intArrayFromList(mapping.getAncestor(helix));
					Arrays.sort(basesInHelixArray);
					double l = computeLengthIncreaseFactor(basesInHelixArray, helix);
					lengthIncreaseFactor.put(helix, l);
					if (l > maxLengthIncreaseFactor) {
						maxLengthIncreaseFactor = l;
						maxIncreaseHelix = helix;
					}
				}
			}
			
			// debug
			System.out.println("Max helix length increase factor = " + maxLengthIncreaseFactor + " reached with helix " + maxIncreaseHelix);;
			
			globalIncreaseFactor = Math.max(1, maxLengthIncreaseFactor);
			
		} else if (helixLengthAdjustmentMethod == DrawRNATemplateMethod.HELIXTRANSLATE) {
			try {
				// Now we need to propagate this helices translations
				Tree<RNANodeValueTemplate> templateAsTree = template.toTree();
				translateVectors = computeHelixTranslations(templateAsTree, mapping);
			
			} catch (ExceptionInvalidRNATemplate e) {
				throw (new RNATemplateDrawingAlgorithmException("ExceptionInvalidRNATemplate: " + e.getMessage()));
			}
		}
		
		// Allocate the coords and centers arrays
		// We create Point2D.Double objects in it but the algorithms
		// we use may choose to create new Point2D.Double objects or to
		// modify those created here.
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()];
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()];
		for (int i = 0; i < _listeBases.size(); i++) {
			coords[i] = new Point2D.Double(0, 0);
			centers[i] = new Point2D.Double(0, 0);
		}
		
		boolean computeCoords = true;
		while (computeCoords) {
			computeCoords = false;
			// Compute coords and centers
			Set<RNATemplateHelix> alreadyDrawnHelixes = new HashSet<RNATemplateHelix>();
			RNATemplateHelix lastMappedHelix = null;
			EdgeEndPoint howWeGotOutOfLastHelix = null;
			int howWeGotOutOfLastHelixBaseIndex = 0;
			iter = template.rnaIterator();
			RNATemplateElement element = null;
			while (iter.hasNext()) {
				element = iter.next();
				if (element instanceof RNATemplateHelix
						&& mapping.getAncestor(element) != null) {
					// We have a mapping between an helix in the RNA sequence
					// and an helix in the template.
					
					RNATemplateHelix helix = (RNATemplateHelix) element;
					boolean firstTimeWeMeetThisHelix;
					int[] basesInHelixArray = RNATemplateAlign.intArrayFromList(mapping.getAncestor(helix));
					Arrays.sort(basesInHelixArray);
					
					// Draw this helix if it has not already been done
					if (!alreadyDrawnHelixes.contains(helix)) {
						firstTimeWeMeetThisHelix = true;
						drawHelixLikeTemplateHelix(basesInHelixArray, helix, coords, centers, globalIncreaseFactor, translateVectors);
						alreadyDrawnHelixes.add(helix);
					} else {
						firstTimeWeMeetThisHelix = false;
					}
					
					EdgeEndPoint howWeGetInCurrentHelix;
					if (firstTimeWeMeetThisHelix) {
						if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
							howWeGetInCurrentHelix = helix.getIn1();
						} else {
							howWeGetInCurrentHelix = helix.getIn2();
						}
					} else {
						if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
							howWeGetInCurrentHelix = helix.getIn2();
						} else {
							howWeGetInCurrentHelix = helix.getIn1();
						}
					}
					
					Point2D.Double P0 = new Point2D.Double();
					Point2D.Double P3 = new Point2D.Double();
					
					if (lastMappedHelix != null) {
						// Now draw the RNA sequence (possibly containing helixes)
						// between the last template drawn helix and this one.
	
						if (lastMappedHelix == helix) {
							// Last helix is the same as the current one so
							// nothing matched (or at best a single
							// non-paired sequence) so we will just
							// use the Radiate algorithm
							
							Point2D.Double helixVector = new Point2D.Double();
							computeHelixEndPointDirections(howWeGotOutOfLastHelix, helixVector, new Point2D.Double());
							
							double angle = angleFromVector(helixVector);
							int b1 = basesInHelixArray[basesInHelixArray.length/2 - 1];
							P0.setLocation(coords[b1]);
							int b2 = basesInHelixArray[basesInHelixArray.length/2];
							P3.setLocation(coords[b2]);
							Point2D.Double loopCenter = new Point2D.Double((P0.x + P3.x)/2, (P0.y + P3.y)/2);
							drawLoop(b1,
									 b2,
									 loopCenter.x,
									 loopCenter.y,
									 angle,
									 coords,
									 centers);
							// If the helix is flipped, we need to compute the symmetric
							// of the whole loop.
							if (helix.isFlipped()) {
								Point2D.Double[] points1 = new Point2D.Double[b2-b1+1];
								Point2D.Double[] points2 = new Point2D.Double[b2-b1+1];
								for (int b=b1; b<=b2; b++) {
									// This is an assignment by reference.
									points1[b-b1] = coords[b];
									points2[b-b1] = centers[b];
								}
								symmetric(loopCenter, helixVector, points1);
								symmetric(loopCenter, helixVector, points2);
							}
						} else {
							// No helices matched between the last helix and
							// the current one, so we draw what is between
							// using the radiate algorithm but on the Bezier curve.
							
							int b1 = howWeGotOutOfLastHelixBaseIndex;
							int b2 = firstTimeWeMeetThisHelix ? basesInHelixArray[0] : basesInHelixArray[basesInHelixArray.length/2];
							P0.setLocation(coords[b1]);
							P3.setLocation(coords[b2]);
							
							Point2D.Double P1, P2;
							
							if (howWeGotOutOfLastHelix.getOtherElement() instanceof RNATemplateUnpairedSequence
									&& howWeGetInCurrentHelix.getOtherElement() instanceof RNATemplateUnpairedSequence) {
								// We will draw the bases on a Bezier curve
								P1 = new Point2D.Double();
								computeBezierTangentVectorTarget(howWeGotOutOfLastHelix, P0, P1);
								
								P2 = new Point2D.Double();
								computeBezierTangentVectorTarget(howWeGetInCurrentHelix, P3, P2);
							} else {
								// We will draw the bases on a straight line between P0 and P3
								P1 = null;
								P2 = null;
							}
							
							drawAlongCurve(b1+1, b2-1, P0, P1, P2, P3, coords, centers);
						}
						
					} else if (basesInHelixArray[0] > 0) {
						// Here we draw what is before the first mapped helix.
	
						RNATemplateUnpairedSequence templateSequence;
						// Try to find our template sequence as the mapped element of base 0
						RNATemplateElement templateSequenceCandidate = mapping.getPartner(0);
						if (templateSequenceCandidate != null
								&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
							templateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
						} else {
							// Try other idea: first template element if it is a sequence
							templateSequenceCandidate = template.getFirst();
							if (templateSequenceCandidate != null
									&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
								templateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
							} else {
								// We don't know where to start
								templateSequence = null;
							}
						}
						
						int b1 = 0;
						int b2 = firstTimeWeMeetThisHelix ? basesInHelixArray[0] : basesInHelixArray[basesInHelixArray.length/2];
						P3.setLocation(coords[b2]);
						
						if (templateSequence != null) {
							coords[0].setLocation(templateSequence.getVertex5());
							coords[0].x *= globalIncreaseFactor;
							coords[0].y *= globalIncreaseFactor;
						} else {
							// We start "near".
							// This is quite arbitrary, the user will need to manually place this bases.
							coords[0].setLocation(coords[b2].x, coords[b2].y + 100);
						}
						P0.setLocation(coords[0]);
						
						Point2D.Double P1, P2;
						
						if (howWeGetInCurrentHelix.getOtherElement() instanceof RNATemplateUnpairedSequence
								&& templateSequence != null) {
							// We will draw the bases on a Bezier curve
							P1 = new Point2D.Double();
							computeBezierTangentVectorTarget(templateSequence.getIn(), P0, P1);
							
							P2 = new Point2D.Double();
							computeBezierTangentVectorTarget(howWeGetInCurrentHelix, P3, P2);
						} else {
							// We will draw the bases on a straight line between P0 and P3
							P1 = null;
							P2 = null;
						}
						
						drawAlongCurve(b1, b2-1, P0, P1, P2, P3, coords, centers);
					}
					
					lastMappedHelix = helix;
					howWeGotOutOfLastHelix = howWeGetInCurrentHelix.getNextEndPoint();
					if (firstTimeWeMeetThisHelix) {
						howWeGotOutOfLastHelixBaseIndex = basesInHelixArray[basesInHelixArray.length/2-1];
					} else {
						howWeGotOutOfLastHelixBaseIndex = basesInHelixArray[basesInHelixArray.length-1];
					}
				}
			} // end template iteration
			
			
			// Now we need to draw what is after the last mapped helix.
			if (howWeGotOutOfLastHelixBaseIndex < coords.length-1
					&& element != null
					&& coords.length > 1) {
				
				RNATemplateUnpairedSequence beginTemplateSequence = null;
				if (lastMappedHelix == null) {
					// No helix at all matched between the template and RNA!
					// So the sequence we want to draw is the full RNA.
					
					// Try to find our template sequence as the mapped element of base 0
					RNATemplateElement templateSequenceCandidate = mapping.getPartner(0);
					if (templateSequenceCandidate != null
							&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
						beginTemplateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
					} else {
						// Try other idea: first template element if it is a sequence
						templateSequenceCandidate = template.getFirst();
						if (templateSequenceCandidate != null
								&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
							beginTemplateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
						} else {
							// We don't know where to start
							beginTemplateSequence = null;
						}
					}
					
					if (beginTemplateSequence != null) {
						coords[0].setLocation(beginTemplateSequence.getVertex5());
						coords[0].x *= globalIncreaseFactor;
						coords[0].y *= globalIncreaseFactor;
					}
					
				}
				
				RNATemplateUnpairedSequence endTemplateSequence;
				// Try to find our template sequence as the mapped element of last base
				RNATemplateElement templateSequenceCandidate = mapping.getPartner(coords.length-1);
				if (templateSequenceCandidate != null
						&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
					endTemplateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
				} else {
					// Try other idea: last template element if it is a sequence
					templateSequenceCandidate = element;
					if (templateSequenceCandidate != null
							&& templateSequenceCandidate instanceof RNATemplateUnpairedSequence) {
						endTemplateSequence = (RNATemplateUnpairedSequence) templateSequenceCandidate;
					} else {
						// We don't know where to end
						endTemplateSequence = null;
					}
				}
				
				int b1 = howWeGotOutOfLastHelixBaseIndex;
				int b2 = coords.length - 1;
				
				if (endTemplateSequence != null) {
					coords[b2].setLocation(endTemplateSequence.getVertex3());
				} else {
					// We end "near".
					// This is quite arbitrary, the user will need to manually place this bases.
					coords[b2].setLocation(coords[b1].x, coords[b1].y + 100);
				}
				coords[coords.length-1].x *= globalIncreaseFactor;
				coords[coords.length-1].y *= globalIncreaseFactor;
				
				if (lastMappedHelix == null
						&& beginTemplateSequence == null
						&& endTemplateSequence == null) {
					// We are extremely unlucky, as we haven't mapped
					// any helix not any half-unpaired sequence from the template!
					
					// Set last base position at something different
					// than (0,0) to prevent all bases to be drawn at (0,0).
					coords[coords.length-1].setLocation(1000, 1000);
				}
				
				Point2D.Double P0 = new Point2D.Double();
				Point2D.Double P3 = new Point2D.Double();
				
				P0.setLocation(coords[b1]);
				P3.setLocation(coords[b2]);
				
				Point2D.Double P1, P2;
				
				if (howWeGotOutOfLastHelix != null
						&& howWeGotOutOfLastHelix.getOtherElement() instanceof RNATemplateUnpairedSequence
						&& endTemplateSequence != null) {
					// We will draw the bases on a Bezier curve
					P1 = new Point2D.Double();
					computeBezierTangentVectorTarget(howWeGotOutOfLastHelix, P0, P1);
					
					P2 = new Point2D.Double();
					computeBezierTangentVectorTarget(endTemplateSequence.getOut(), P3, P2);
				} else if (lastMappedHelix == null
						&& beginTemplateSequence != null
						&& endTemplateSequence != null) {
					// We will draw the bases on a Bezier curve
					P1 = new Point2D.Double();
					computeBezierTangentVectorTarget(beginTemplateSequence.getIn(), P0, P1);
					
					P2 = new Point2D.Double();
					computeBezierTangentVectorTarget(endTemplateSequence.getOut(), P3, P2);
				} else {
					// We will draw the bases on a straight line between P0 and P3
					P1 = null;
					P2 = null;
				}
				
				drawAlongCurve((lastMappedHelix != null ? b1+1 : b1), b2, P0, P1, P2, P3, coords, centers);
			
			}
			
			
			if (helixLengthAdjustmentMethod == DrawRNATemplateMethod.NOINTERSECT && coords.length > 3) {
				// Are we happy with this value of globalIncreaseFactor?
				Line2D.Double[] lines = new Line2D.Double[coords.length-1];
				for (int i=0; i<coords.length-1; i++) {
					lines[i] = new Line2D.Double(coords[i], coords[i+1]);
				}
				int intersectLines = 0;
				for (int i=0; i<lines.length; i++) {
					for (int j=i+2; j<lines.length; j++) {
						if (lines[i].intersectsLine(lines[j])) {
							intersectLines++;
						}
					}
				}
				// If no intersection we keep this globalIncreaseFactor value
				if (intersectLines > 0) {
					// Don't increase more than a maximum value
					if (globalIncreaseFactor < 3) {
						globalIncreaseFactor += 0.1;
						System.out.println("globalIncreaseFactor increased to " + globalIncreaseFactor);
						// Compute the drawing again
						computeCoords = true;
					}
				}
			}
			
		}
		
		// debug
		if (helixLengthAdjustmentMethod == DrawRNATemplateMethod.MAXSCALINGFACTOR
				|| helixLengthAdjustmentMethod == DrawRNATemplateMethod.NOINTERSECT) {
			System.out.println("globalIncreaseFactor = " + globalIncreaseFactor);
		}
		
		// Now we actually move the bases, according to arrays coords and centers
		// and taking in account the space between bases parameter.
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(coords[i].x * _spaceBetweenBases,
							coords[i].y * _spaceBetweenBases));
			_listeBases.get(i).setCenter(
					new Point2D.Double(centers[i].x * _spaceBetweenBases,
							centers[i].y * _spaceBetweenBases));
		}
		return mapping;
	}
	
	



	private static  double objFun(int n1, int n2, double r, double bpdist, double multidist) {
		return (((double) n1) * 2.0
				* Math.asin(((double) bpdist) / (2.0 * r))
				+ ((double) n2) * 2.0
				* Math.asin(((double) multidist) / (2.0 * r)) - (2.0 * Math.PI));
	}

	public double determineRadius(int nbHel, int nbUnpaired, double startRadius) {
	  return determineRadius(nbHel,nbUnpaired,startRadius,BASE_PAIR_DISTANCE,MULTILOOP_DISTANCE);
	}

	
	public static double determineRadius(int nbHel, int nbUnpaired, double startRadius, double bpdist, double multidist) {
		double xmin = bpdist / 2.0;
		double xmax = 3.0 * multidist + 1;
		double x = (xmin + xmax) / 2.0;
		double y = 10000.0;
		double ymin = -1000.0;
		double ymax = 1000.0;
		int numIt = 0;
		double precision = 0.00001;
		while ((Math.abs(y) > precision) && (numIt < 10000)) {
			x = (xmin + xmax) / 2.0;
			y = objFun(nbHel, nbUnpaired, x, bpdist, multidist);
			ymin = objFun(nbHel, nbUnpaired, xmax, bpdist, multidist);
			ymax = objFun(nbHel, nbUnpaired, xmin, bpdist, multidist);
			if (ymin > 0.0) {
				xmax = xmax + (xmax - xmin);
			} else if ((y <= 0.0) && (ymax > 0.0)) {
				xmax = x;
			} else if ((y >= 0.0) && (ymin < 0.0)) {
				xmin = x;
			} else if (ymax < 0.0) {
				xmin = Math.max(xmin - (x - xmin), Math.max(
						bpdist / 2.0, multidist / 2.0));
				xmax = x;
			}
			numIt++;
		}
		return x;
	}

	public void drawRNA(VARNAConfig conf) throws ExceptionNAViewAlgorithm {
		drawRNA(RNA.DEFAULT_DRAW_MODE,conf);
	}
	public void drawRNA(int mode,VARNAConfig conf) throws ExceptionNAViewAlgorithm {
		_drawMode = mode;
		switch (get_drawMode()) {
		case RNA.DRAW_MODE_RADIATE:
			drawRNARadiate(conf);
			break;
		case RNA.DRAW_MODE_LINEAR:
			drawRNALine();
			break;
		case RNA.DRAW_MODE_CIRCULAR:
			drawRNACircle();
			break;
		case RNA.DRAW_MODE_NAVIEW:
			drawRNANAView();
			break;
		case RNA.DRAW_MODE_VARNA_VIEW:
			drawRNAVARNAView();
			break;
		case RNA.DRAW_MODE_MOTIFVIEW:
			drawMOTIFView();
			break;
		default:
			break;
		}

	}
	

	public int getDrawMode() {
		return _drawMode;
	}

	
	
	private void drawLoop(int i, int j, double x, double y, double dirAngle,
			Point2D.Double[] coords, Point2D.Double[] centers) {
		if (i > j) {
			return;
		}

		// BasePaired
		if (_listeBases.get(i).getElementStructure() == j) {
			double normalAngle = Math.PI / 2.0;
			centers[i] = new Point2D.Double(x, y);
			centers[j] = new Point2D.Double(x, y);
			coords[i].x = (x + BASE_PAIR_DISTANCE
					* Math.cos(dirAngle - normalAngle) / 2.0);
			coords[i].y = (y + BASE_PAIR_DISTANCE
					* Math.sin(dirAngle - normalAngle) / 2.0);
			coords[j].x = (x + BASE_PAIR_DISTANCE
					* Math.cos(dirAngle + normalAngle) / 2.0);
			coords[j].y = (y + BASE_PAIR_DISTANCE
					* Math.sin(dirAngle + normalAngle) / 2.0);
			drawLoop(i + 1, j - 1, x + LOOP_DISTANCE * Math.cos(dirAngle), y
					+ LOOP_DISTANCE * Math.sin(dirAngle), dirAngle, coords,
					centers);
		} else {
			int k = i;
			Vector<Integer> basesMultiLoop = new Vector<Integer>();
			Vector<Integer> helices = new Vector<Integer>();
			int l;
			while (k <= j) {
				l = _listeBases.get(k).getElementStructure();
				if (l > k) {
					basesMultiLoop.add(new Integer(k));
					basesMultiLoop.add(new Integer(l));
					helices.add(new Integer(k));
					k = l + 1;
				} else {
					basesMultiLoop.add(new Integer(k));
					k++;
				}
			}
			int mlSize = basesMultiLoop.size() + 2;
			int numHelices = helices.size() + 1;
			double totalLength = MULTILOOP_DISTANCE * (mlSize - numHelices)
					+ BASE_PAIR_DISTANCE * numHelices;
			double multiLoopRadius;
			double angleIncrementML;
			double angleIncrementBP;
			if (mlSize > 3) {
				multiLoopRadius = determineRadius(numHelices, mlSize
						- numHelices, (totalLength) / (2.0 * Math.PI),
						BASE_PAIR_DISTANCE,
						MULTILOOP_DISTANCE
						);
				angleIncrementML = -2.0
						* Math.asin(((float) MULTILOOP_DISTANCE)
								/ (2.0 * multiLoopRadius));
				angleIncrementBP = -2.0
						* Math.asin(((float) BASE_PAIR_DISTANCE)
								/ (2.0 * multiLoopRadius));
			} else {
				multiLoopRadius = 35.0;
				angleIncrementBP = -2.0
						* Math.asin(((float) BASE_PAIR_DISTANCE)
								/ (2.0 * multiLoopRadius));
				angleIncrementML = (-2.0 * Math.PI - angleIncrementBP) / 2.0;
			}
			// System.out.println("MLr:"+multiLoopRadius+" iBP:"+angleIncrementBP+" iML:"+angleIncrementML);

			double centerDist = Math.sqrt(Math.max(Math.pow(multiLoopRadius, 2)
					- Math.pow(BASE_PAIR_DISTANCE / 2.0, 2), 0.0))
					- LOOP_DISTANCE;
			Point2D.Double mlCenter = new Point2D.Double(
					(x + (centerDist * Math.cos(dirAngle))),
					(y + (centerDist * Math.sin(dirAngle))));

			// Base directing angle for (multi|hairpin) loop, from the center's
			// perspective
			double baseAngle = dirAngle
			// U-turn
					+ Math.PI
					// Account for already drawn supporting base-pair
					+ 0.5 * angleIncrementBP
					// Base cannot be paired twice, so next base is at
					// "unpaired base distance"
					+ 1.0 * angleIncrementML;
			double[] angles = new double[_listeBases.size()];
			int n1 = 1;
			int n2 = 1;
			for (k = basesMultiLoop.size() - 1; k >= 0; k--) {
				l = basesMultiLoop.get(k).intValue();
				centers[l] = mlCenter;
				angles[l] = baseAngle;
				coords[l].x = mlCenter.x + multiLoopRadius
						* Math.cos(baseAngle);
				coords[l].y = mlCenter.y + multiLoopRadius
						* Math.sin(baseAngle);
				if ((_listeBases.get(l).getElementStructure() < l)
						&& (_listeBases.get(l).getElementStructure() != -1)) {
					baseAngle += angleIncrementBP;
					n1++;
				} else {
					baseAngle += angleIncrementML;
					n2++;
				}
			}
			// System.out.println("n1:"+n1+" n2:"+n2);
			double newAngle;
			int m, n;
			for (k = 0; k < helices.size(); k++) {
				m = helices.get(k).intValue();
				n = _listeBases.get(m).getElementStructure();
				newAngle = (angles[m] + angles[n]) / 2.0;
				drawLoop(m + 1, n - 1, (LOOP_DISTANCE * Math.cos(newAngle))
						+ (coords[m].x + coords[n].x) / 2.0,
						(LOOP_DISTANCE * Math.sin(newAngle))
								+ (coords[m].y + coords[n].y) / 2.0, newAngle,
						coords, centers);
			}
		}
	}


	public void drawRNARadiate(VARNAConfig conf) {
		drawRNARadiate(-1.0,conf._flatExteriorLoop);
	}
	
	public void drawRNARadiate(double dirAngle, boolean flatExteriorLoop) {
		_drawn = true;
		_drawMode = DRAW_MODE_RADIATE;
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()];
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()];
		for (int i = 0; i < _listeBases.size(); i++) {
			coords[i] = new Point2D.Double(0, 0);
			centers[i] = new Point2D.Double(0, 0);
		}
		if (flatExteriorLoop)
		{
		  dirAngle += 1.0 - Math.PI/2.0;
		  int i=0;
		  double x = 0.0;
		  double y = 0.0;
		  double vx = -Math.sin(dirAngle);
		  double vy = Math.cos(dirAngle);
		  while(i<_listeBases.size())
		  {
			  coords[i].x = x; 
			  coords[i].y = y;
			  centers[i].x = x+BASE_PAIR_DISTANCE*vy;
			  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			  int j = _listeBases.get(i).getElementStructure();
			  if (j>i)
			  {
				  drawLoop(i, j, x+(BASE_PAIR_DISTANCE*vx/2.0), y+(BASE_PAIR_DISTANCE*vy/2.0), dirAngle, coords, centers);
				  centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
				  i = j;
				  x += BASE_PAIR_DISTANCE*vx;
				  y += BASE_PAIR_DISTANCE*vy;
				  centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			  }
			  x += MULTILOOP_DISTANCE*vx;
			  y += MULTILOOP_DISTANCE*vy;
			  i += 1;
		  }
		}
		else
		{
		  drawLoop(0, _listeBases.size() - 1, 0, 0, dirAngle, coords, centers);
		}
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(coords[i].x * _spaceBetweenBases,
							coords[i].y * _spaceBetweenBases));
			_listeBases.get(i).setCenter(
					new Point2D.Double(centers[i].x * _spaceBetweenBases,
							centers[i].y * _spaceBetweenBases));
		}

		// TODO
		// change les centres des bases de la premiere helice vers la boucle la
		// plus proche
	}

	public void drawRNANAView() throws ExceptionNAViewAlgorithm {
		_drawMode = DRAW_MODE_NAVIEW;
		_drawn = true;

		ArrayList<Double> X = new ArrayList<Double>(_listeBases.size());
		ArrayList<Double> Y = new ArrayList<Double>(_listeBases.size());
		ArrayList<Short> pair_table = new ArrayList<Short>(_listeBases.size());

		for (int i = 0; i < _listeBases.size(); i++) {
			pair_table.add(Short.valueOf(String.valueOf(_listeBases.get(i)
					.getElementStructure())));
		}
		NAView naView = new NAView();
		naView.naview_xy_coordinates(pair_table, X, Y);

		// Updating individual base positions
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(X.get(i) * 2.5 * _spaceBetweenBases, Y
							.get(i)
							* 2.5 * _spaceBetweenBases));
		}

		// Updating centers
		for (int i = 0; i < _listeBases.size(); i++) {
			int indicePartner = _listeBases.get(i).getElementStructure();
			if (indicePartner != -1) {
				Point2D.Double base = _listeBases.get(i).getCoords();
				Point2D.Double partner = _listeBases.get(indicePartner)
						.getCoords();
				_listeBases.get(i).setCenter(
						new Point2D.Double((base.x + partner.x) / 2.0,
								(base.y + partner.y) / 2.0));
			} 
			else {
				Vector<Integer> loop = getLoopBases(i);
				double tmpx = 0.0;
				double tmpy = 0.0;
				for (int j = 0; j < loop.size(); j++) {
					int partner = loop.elementAt(j);
					Point2D.Double loopmember = _listeBases.get(partner)
							.getCoords();
					tmpx += loopmember.x;
					tmpy += loopmember.y;
				}
				_listeBases.get(i).setCenter(
						new Point2D.Double(tmpx / loop.size(), tmpy
								/ loop.size()));
			}
		}
	}

	public void drawMOTIFView() {
		_drawn = true;
		_drawMode = DRAW_MODE_MOTIFVIEW;	
		int spaceBetweenStrand =0;
		Motif motif = new Motif(this,get_listeBases());
		motif.listStrand();
		for (int i = 0; i < motif.getListStrand().sizeStruct(); i++ ){
			for (int j = 0; j < motif.getListStrand().getStrand(i).sizeStrand(); j++ ){
				int indice = motif.getListStrand().getStrand(i).getMB(j).getIndex();
				get_listeBases().get(indice).setCoords(
						new Point2D.Double(0,0));
				get_listeBases().get(indice).setCenter(
						new Point2D.Double(0, 0));

			}
		}
		//Recherche du brin central
		int centralStrand = motif.getCentralStrand();
		
		//Cas o l'on a un motif en toile
		if(centralStrand!=-1){
			//On positionne le brin central
			motif.positionneSpecificStrand(centralStrand, spaceBetweenStrand);

			//On place les autres brins par rapport a ce brin central
			motif.orderStrands(centralStrand);
		}
		
		else {	
			centralStrand = 0;
			motif.positionneStrand();
			motif.ajusteStrand();
		}
		motif.reajustement();
		motif.deviationBasePair();	
		motif.setCenterMotif();
	}
	
		
	public ArrayList<ModeleBase> getAllPartners(int indice) {		
		ArrayList<ModeleBase> result = new ArrayList<ModeleBase>();
		ModeleBase me = this.getBaseAt(indice);
		int i = me.getElementStructure();
		if (i!= -1)
		{ result.add(getBaseAt(i)); }
		ArrayList<ModeleBP> msbps = getAuxBPs(indice);
		for (ModeleBP m : msbps)
		{ result.add(m.getPartner(me)); }		
		return result;
	}
	
	
	public int get_drawMode() {
		return _drawMode;
	}

	public void setDrawMode(int drawMode) {
		_drawMode = drawMode;
	}

	public void setRNA(String seq, String str)
		throws ExceptionFileFormatOrSyntax, ExceptionUnmatchedClosingParentheses 
	{
		setRNA(RNA.explodeSequence(seq), str);
	}

	
	public void setRNA(String[] seq, int[] str)
		throws ExceptionFileFormatOrSyntax 
	{
		setRNA(seq, str, 1);
	}

	public void setRNA(List<String> seq, int[] str)
	throws ExceptionFileFormatOrSyntax 
	{
		setRNA(seq.toArray(new String[seq.size()]), str, 1);
	}
	
	public void setRNA(List<String> seq, int[] str, int baseIndex) throws ExceptionFileFormatOrSyntax
	{
		setRNA(seq.toArray(new String[seq.size()]), str, baseIndex);
	}

	
	public void setRNA(String[] seq, int[] str, int baseIndex)
			throws ExceptionFileFormatOrSyntax {
		clearAnnotations();
		_listeBases = new ArrayList<ModeleBase>();
		if (seq.length != str.length) {
			warningEmition("Sequence length " + seq.length
					+ " differs from that of secondary structure " + str.length
					+ ". \nAdapting sequence length ...");
			if (seq.length < str.length) {
				String [] nseq = new String[str.length];
				for(int i=0;i<seq.length;i++)
				{
					nseq[i] = seq[i];
				}
				for(int i=seq.length;i<nseq.length;i++)
				{
					nseq[i] = "";
				}
				seq = nseq;
			} else {
				String[] seqTmp = new String[str.length];
				for (int i =0;i<str.length;i++)
				{
					seqTmp[i] = seq[i]; 
				}
				seq = seqTmp; 
			}
		}
		for (int i = 0; i < str.length; i++) {
			_listeBases.add(new ModeleBaseNucleotide(seq[i], i, baseIndex+i));
		}
		applyStruct(str);
	}

	/**
	 * Sets the RNA to be drawn. Uses when comparison mode is on. Will draw the
	 * super-structure passed in parameters and apply specials styles to the
	 * bases owning by each RNA alignment and both.
	 * 
	 * @param seq
	 *            - The sequence of the super-structure This sequence shall be
	 *            designed like this:
	 *            <code>firstRNA1stBaseSecondRNA1stBaseFirstRNA2ndBaseSecondRNA2ndBase [...]</code>
	 * <br>
	 *            <b>Example:</b> <code>AAC-GUAGA--UGG</code>
	 * @param struct
	 *            - The super-structure
	 * @param basesOwn
	 *            - The RNA owning bases array (each index will be:0 when common
	 *            base, 1 when first RNA alignment base, 2 when second RNA
	 *            alignment base)
	 * @throws ExceptionUnmatchedClosingParentheses
	 * @throws ExceptionFileFormatOrSyntax 
	 */
	public void setRNA(String seq, String struct, ArrayList<Integer> basesOwn)
			throws ExceptionUnmatchedClosingParentheses, ExceptionFileFormatOrSyntax {
		clearAnnotations();
		_listeBases = new ArrayList<ModeleBase>();
		// On "parse" la structure (repérage des points, tiret et couples
		// parentheses ouvrante/fermante)
		int[] array_struct = parseStruct(struct);
		int size = struct.length();
		int j = 0;
		for (int i = 0; i < size; i++) {
			ModeleBase mb;
			if (seq.charAt(j)!=seq.charAt(j + 1))
			{
				ModeleBasesComparison mbc = new ModeleBasesComparison(seq.charAt(j), seq.charAt(j + 1), i);
				mbc.set_appartenance(basesOwn.get(i));
				mbc.setBaseNumber(i+1);
				mb = mbc;
			}
			else
			{
			    mb = new ModeleBaseNucleotide(""+seq.charAt(j),i,i+1);
			    
			}
    		_listeBases.add(mb);
			j+=2;
		}
		for (int i = 0; i < size; i++) {
			if (array_struct[i]!=-1)
			{this.addBP(i, array_struct[i]);}
			
			j += 2;
		}
	}

	public void setRNA(List<String> seq, String dbnStr)
			throws ExceptionUnmatchedClosingParentheses,
			ExceptionFileFormatOrSyntax {
		clearAnnotations();
		String parDBN = dbnStr.replace('(', '(').replace(')', ')').replace('[', ':').replace(']', ':').replace('{', ':').replace('}', ':').replace('<', ':').replace('>', ':');
		String braDBN = dbnStr.replace('(', ':').replace(')', ':').replace('[', '(').replace(']', ')').replace('{', ':').replace('}', ':').replace('<', ':').replace('>', ':');
		String accDBN = dbnStr.replace('(', ':').replace(')', ':').replace('[', ':').replace(']', ':').replace('{', '(').replace('}', ')').replace('<', ':').replace('>', ':');
		String cheDBN = dbnStr.replace('(', ':').replace(')', ':').replace('[', ':').replace(']', ':').replace('{', ':').replace('}', ':').replace('<', '(').replace('>', ')');
	    int[] parStr = parseStruct(parDBN);
		int[] braStr = parseStruct(braDBN);
		int[] accStr = parseStruct(accDBN);
		int[] cheStr = parseStruct(cheDBN);
		int[] finStr = new int[parStr.length];
		for(int i=0;i<parStr.length;i++)
		{ finStr[i] = -1; } 
		
		for(int i=0;i<parStr.length;i++)
		{
			if (parStr[i]>i)
			{
				finStr[i] = parStr[i];
				finStr[finStr[i]] = i;
			}
			else if (braStr[i]>i)
			{
				if ((parStr[i]==-1)&&(parStr[braStr[i]]==-1))
				{
					finStr[i] = braStr[i];
					finStr[finStr[i]] = i;				
				}
			}
			else if (accStr[i]>i)
			{
				if ((parStr[i]==-1)&&(parStr[accStr[i]]==-1)&&(braStr[i]==-1)&&(braStr[accStr[i]]==-1))
				{
					finStr[i] = accStr[i];
					finStr[finStr[i]] = i;				
				}
			}
			else if (cheStr[i]>i)
			{
				if ((parStr[i]==-1)&&(parStr[cheStr[i]]==-1)&&(braStr[i]==-1)&&(braStr[cheStr[i]]==-1)&&(accStr[i]==-1)&&(accStr[cheStr[i]]==-1))
				{
					finStr[i] = cheStr[i];
					finStr[cheStr[i]] = i;				
				}
			}
		}
		setRNA(seq, finStr);
	}

	public static ArrayList<String> explodeSequence(String seq)
	{
		ArrayList<String> analyzedSeq = new ArrayList<String>();
		int i=0;
		while(i<seq.length())
		{
			if (seq.charAt(i)=='{')
			{
				boolean found = false;
				String buf = "";
				i++;
				while(!found & (i<seq.length()))
				{
					if (seq.charAt(i)!='}')
					{ 
						buf += seq.charAt(i);
						i++;
					}
					else
					{ found = true;}
				}
				analyzedSeq.add(buf);
			}
			else
			{
				analyzedSeq.add(""+seq.charAt(i));
			}
			i++;
		}
		return analyzedSeq;
	}
	
	public int[] parseStruct(String str)
			throws ExceptionUnmatchedClosingParentheses, ExceptionFileFormatOrSyntax {
		int[] result = new int[str.length()];
		int unexpectedChar = -1;
		Stack<Integer> p = new Stack<Integer>();
		for (int i = 0; i < str.length(); i++) {
			char c = str.charAt(i);
			if (c == '(') {
				p.push(new Integer(i));
			} else if (c == '.' || c == '-' || c == ':') {
				result[i] = -1;
			} else if (c == ')') {
				if (p.size() == 0) {
					throw new ExceptionUnmatchedClosingParentheses(i + 1);
				}
				int j = p.pop().intValue();
				result[i] = j;
				result[j] = i;
			} else {
				if (unexpectedChar == -1)
					unexpectedChar = i;
				break;
			}
		}

		if (unexpectedChar != -1) {
			//warningEmition("Unexpected Character at index:" + unexpectedChar);
		}

		if (p.size() != 0) {
			throw new ExceptionUnmatchedClosingParentheses(
					p.pop().intValue() + 1);
		}

		return result;
	}

	public Point getHelixInterval(int index) {
		if ((index < 0) || (index >= _listeBases.size())) {
			return new Point(index, index);
		}
		int j = _listeBases.get(index).getElementStructure();
		if (j != -1) {
			int minH = index;
			int maxH = index;
			if (j > index) {
				maxH = j;
			} else {
				minH = j;
			}
			boolean over = false;
			while (!over) {
				if ((minH < 0) || (maxH >= _listeBases.size())) {
					over = true;
				} else {
					if (_listeBases.get(minH).getElementStructure() == maxH) {
						minH--;
						maxH++;
					} else {
						over = true;
					}
				}
			}
			minH++;
			maxH--;
			return new Point(minH, maxH);
		}
		return new Point(0, 0);
	}

	public ArrayList<Integer> getHelix(int index) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		if ((index < 0) || (index >= _listeBases.size())) {
			return result;
		}
		Point p = getHelixInterval(index);
		for (int i =p.x;i<=p.y;i++)
		{
			result.add(i);
			result.add(this._listeBases.get(i).getElementStructure());
		}
		return result;
	}

	public Point getMultiLoop(int index) {
		if ((index < 0) || (index >= _listeBases.size())) {
			return new Point(index, index);
		}
		Point h = getHelixInterval(index);
		int minH = h.x - 1;
		int maxH = h.y + 1;
		boolean over = false;
		while (!over) {
			if (minH < 0) {
				over = true;
				minH = 0;
			} else {
				if (_listeBases.get(minH).getElementStructure() == -1) {
					minH--;
				} else if (_listeBases.get(minH).getElementStructure() < minH) {
					minH = _listeBases.get(minH).getElementStructure() - 1;
				} else {
					over = true;
				}
			}
		}
		over = false;
		while (!over) {
			if (maxH > _listeBases.size() - 1) {
				over = true;
				maxH = _listeBases.size() - 1;
			} else {
				if (_listeBases.get(maxH).getElementStructure() == -1) {
					maxH++;
				} else if (_listeBases.get(maxH).getElementStructure() > maxH) {
					maxH = _listeBases.get(maxH).getElementStructure() + 1;
				} else {
					over = true;
				}
			}
		}
		return new Point(minH, maxH);
	}

	public Vector<Integer> getLoopBases(int startIndex) {
		Vector<Integer> result = new Vector<Integer>();

		if ((startIndex < 0) || (startIndex >= _listeBases.size())) {
			return result;
		}
		int index = startIndex;
		result.add(startIndex);
		if (_listeBases.get(index).getElementStructure() <= index) {
			index = (index + 1) % _listeBases.size();
		} else {
			index = _listeBases.get(index).getElementStructure();
			result.add(index);
			index = (index + 1) % _listeBases.size();
		}

		while (index != startIndex) {
			result.add(index);
			if (_listeBases.get(index).getElementStructure() == -1) {
				index = (index + 1) % _listeBases.size();
			} else {
				index = _listeBases.get(index).getElementStructure();
				result.add(index);
				index = (index + 1) % _listeBases.size();
			}
		}
		return result;
	}

	/**
	 * Returns the RNA secondary structure displayed by this panel as a
	 * well-parenthesized word, accordingly to the DBN format
	 * 
	 * @return This panel's secondary structure
	 */
	public String getStructDBN() {
		String result = "";
		for (int i = 0; i < _listeBases.size(); i++) {
			int j = _listeBases.get(i).getElementStructure();
			if (j == -1) {
				result += ".";
			} else if (i > j) {
				result += ")";
			} else {
				result += "(";
			}
		}
		return result;
	}
	
	public ArrayList<int[]> paginateStructure()
	{
		ArrayList<int[]> result = new ArrayList<int[]>();
		
		// Mumbo jumbo to sort the basepair list
		ArrayList<ModeleBP> bps =  this.getAllBPs();
		ModeleBP[] mt = new ModeleBP[bps.size()];
		bps.toArray(mt);
		Arrays.sort(mt, new Comparator<ModeleBP>(){
			public int compare(ModeleBP arg0, ModeleBP arg1) {
				if (arg0.getIndex5()!=arg1.getIndex5())
					return arg0.getIndex5()-arg1.getIndex5();
				else
					return arg0.getIndex3()-arg1.getIndex3();
				
			}});
		bps = new ArrayList<ModeleBP>();
		for(int i=0;i<mt.length;i++)
		{
			bps.add(mt[i]);
		}
		
		
		while (bps.size()!=0)
		{
			ArrayList<ModeleBP> currentBPs = new ArrayList<ModeleBP>();
			Stack<Integer> pile = new Stack<Integer>();
			int[] ss = new int[this.getSize()];
			for (int i=0;i<ss.length;i++)
			{ ss[i] = -1; }
			for(int i=0;i<bps.size();i++)
			{
				ModeleBP bp = bps.get(i);
				boolean ok = true;
				if (!pile.empty())
				{
				  int x = pile.peek();
				  if ((bp.getIndex5()<=x && bp.getIndex3()>=x) || (ss[bp.getIndex5()]!=-1))
				  {
					  ok = false;
				  }
				}
				if (ok)
				{
					ss[bp.getIndex5()] = bp.getIndex3();
					ss[bp.getIndex3()] = bp.getIndex5();
					currentBPs.add(bp);
					pile.add(bp.getIndex3());
				}
				if (!pile.empty() && (i==pile.peek()))
				{
					pile.pop();
				}
			}
			bps.removeAll(currentBPs);
			result.add(ss);
		}
		return result;
	}
	
	
	public String getStructDBN(boolean includeMostPKs) {
		String result = getStructDBN();
		if (includeMostPKs)
		{
			ArrayList<int[]> pages = paginateStructure();
			char[] res = new char[getSize()]; 
			for (int i=0;i<res.length;i++)
			{
				res[i] = '.';
			}
			System.out.println(pages.size());
			char[] open = {'(','[','{','<'};
			char[] close = {')',']','}','>'};
			for (int p=0;p<Math.min(pages.size(),open.length);p++)
			{
				int[] page = pages.get(p);
				for (int i=0;i<res.length;i++)
				{
					if (page[i]!=-1 && page[i]>i && res[i]=='.' && res[page[i]]=='.')
					{
						res[i]=open[p];
						res[page[i]]=close[p];
					}
				}
			}
			result = "";
			for (int i=0;i<res.length;i++)
			{
				result += res[i];
			}

		}
		return result;
		
	}

	public String getStructDBN(int[] str) {
		String result = "";
		for (int i = 0; i < str.length; i++) {
			if (str[i] == -1) {
				result += ".";
			} else if (str[i] > i) {
				result += "(";
			} else {
				result += ")";
			}
		}
		return result;
	}

	/**
	 * Returns the raw nucleotides sequence for the displayed RNA
	 * 
	 * @return The RNA sequence
	 */
	public String getSeq() {
		String result = "";
			for (int i = 0; i < _listeBases.size(); i++) {
				result += ((ModeleBase) _listeBases.get(i)).getContent();
			}
		return result;
	}

	public String getStructBPSEQ() {
		String result = "";
			int[] str = getNonOverlappingStruct();
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1) + " "
						+ ((ModeleBaseNucleotide) _listeBases.get(i)).getContent()
						+ " " + (str[i] + 1)
						+ "\n";
			}
		return result;
	}
	
	public int[] getNonCrossingStruct()
	{
	  int[] result = new int[_listeBases.size()];
      // Adding "planar" base-pairs
	  for (int i = 0; i < _listeBases.size(); i++) {
			result[i] = _listeBases.get(i).getElementStructure();
		}
	  return result;
	}

	public int[] getNonOverlappingStruct()
	{
	  int[] result = getNonCrossingStruct();
	  // Adding additional base pairs when possible (No more than one base-pair per base)
	  for (int i = 0; i < _structureAux.size(); i++) {
			ModeleBP msbp = _structureAux.get(i);
			ModeleBase mb5 = msbp.getPartner5();
			ModeleBase mb3 = msbp.getPartner3();
			int j5 = mb5.getIndex();
			int j3 = mb3.getIndex();
			if ((result[j3]==-1) && (result[j5]==-1))
			{
				result[j3] = j5;
				result[j5] = j3;
			}
		}
	  return result;
	}
	
	
	public String getStructCT() {
		String result = "";
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1)
						+ " "
						+ ((ModeleBase) _listeBases.get(i))
								.getContent() + " " + i + " " + (i + 2) + " "
						+ (_listeBases.get(i).getElementStructure() + 1) + " "
						+ (i + 1) + "\n";
			}
		return result;
	}

	public void saveAsBPSEQ(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("# " + title + "\n");
			f.write(this.getStructBPSEQ() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public void saveAsCT(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("" + _listeBases.size() + " " + title + "\n");
			f.write(this.getStructCT() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public void saveAsDBN(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("> " + title + "\n");
			f.write(getListeBasesToString() + "\n");
			f.write(getStructDBN() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public String getListeBasesToString() {
		String s = new String();
			for (int i = 0; i < _listeBases.size(); i++) {
				s += ((ModeleBaseNucleotide) _listeBases.get(i)).getContent();
		}
		return s;
	}

	
	public void applyBPs(ArrayList<ModeleBP> allbps)
	{
		ArrayList<ModeleBP> planar = new ArrayList<ModeleBP>(); 
		ArrayList<ModeleBP> others = new ArrayList<ModeleBP>(); 
		//System.err.println("Sequence: "+this.getSeq());
		RNAMLParser.planarize(allbps, planar, others, getSize());
		//System.err.println("All:"+allbps);
		//System.err.println("=> Planar: "+planar);
		//System.err.println("=> Others: "+others);
		
		for (ModeleBP mb :planar) {
			addBP(mb.getPartner5().getIndex(), mb.getPartner3().getIndex(), mb);
		}

		for (ModeleBP mb :others) {
			addBPAux(mb.getPartner5().getIndex(), mb.getPartner3().getIndex(), mb);
		}
	}
	

	public void set_listeBases(ArrayList<ModeleBase> _liste) {
		this._listeBases = _liste;
	}

	public void addVARNAListener(InterfaceVARNAListener rl) {
		_listeVARNAListener.add(rl);
	}

	public void warningEmition(String warningMessage) {
		for (int i = 0; i < _listeVARNAListener.size(); i++) {
			_listeVARNAListener.get(i).onWarningEmitted(warningMessage);
		}
	}

	public void applyStyleOnBases(ArrayList<Integer> basesList,
			ModeleStyleBase style) {
		for (int i = 1; i < basesList.size(); i++) {
			_listeBases.get(basesList.get(i)).setStyleBase(style);
		}
	}

	private int[] correctReciprocity(int[] str) {
		int[] result = new int[str.length];
		for (int i = 0; i < str.length; i++) {
			if (str[i] != -1) {
				if (i == str[str[i]]) {
					result[i] = str[i];
				}
				else {
					str[str[i]] = i;
				}
			}
			else
			{
				result[i] = -1;
			}
		}
		return result;
	}

	private void applyStruct(int[] str) throws ExceptionFileFormatOrSyntax {
		str = correctReciprocity(str);

		int[] planarSubset = RNAMLParser.planarize(str);
		_structureAux.clear();

		for (int i = 0; i < planarSubset.length; i++) {
			if (str[i]>i)
			{
				if (planarSubset[i] > i) 
				{  addBP(i,planarSubset[i]);  }
				else if ((planarSubset[i]!=str[i]))
				{  addBPAux(i,str[i]);  }
			}
		}

	}



	public ArrayList<ModeleBase> get_listeBases() {
		return _listeBases;
	}


	public int getSize() {
		return _listeBases.size();
	}

	public ArrayList<Integer> findAll() {
		ArrayList<Integer> listAll = new ArrayList<Integer>();
		for (int i = 0; i < get_listeBases().size(); i++) {
			listAll.add(i);
		}
		return listAll;
	}

	public ArrayList<Integer> findBulge(int index) {
		ArrayList<Integer> listUp = new ArrayList<Integer>();
		if (get_listeBases().get(index).getElementStructure() == -1) {
			int i = index;
			boolean over = false;
			while ((i < get_listeBases().size()) && !over) {
				int j = get_listeBases().get(i).getElementStructure();
				if (j == -1) {
					listUp.add(i);
					i++;
				} else {
					over = true;
				}
			}
			i = index - 1;
			over = false;
			while ((i >= 0) && !over) {
				int j = get_listeBases().get(i).getElementStructure();
				if (j == -1) {
					listUp.add(i);
					i--;
				} else {
					over = true;
				}
			}
		}
		return listUp;
	}

	public ArrayList<Integer> findStem(int index) {
		ArrayList<Integer> listUp = new ArrayList<Integer>();
		int i = index;
		do {
			listUp.add(i);
			int j = get_listeBases().get(i).getElementStructure();
			if (j == -1) {
				i = (i + 1) % getSize();
			} else {
				if ((j < i) && (index <= i) && (j <= index)) {
					i = j;
				} else {
					i = (i + 1) % getSize();
				}
			}
		} while (i != index);
		return listUp;
	}

	public int getHelixCountOnLoop(int indice) {
		int cptHelice = 0;
		if (indice < 0 || indice >= get_listeBases().size())
			return cptHelice;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if ((j != -1) && (j < i)) {
			i = j + 1;
			indice = i;
		}
		do {
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
				cptHelice++;
			} else {
				i = (i + 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return cptHelice;
	}

	public ArrayList<Integer> findLoop(int indice) {
		return findLoopForward(indice);
	}

	public ArrayList<Integer> findLoopForward(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		if (indice < 0 || indice >= get_listeBases().size())
			return base;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if (j != -1) {
			i = Math.min(i, j) + 1;
			indice = i;
		}
		do {
			base.add(i);
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
			} else {
				i = (i + 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return base;
	}

	public ArrayList<Integer> findPair(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		int j = get_listeBases().get(indice).getElementStructure();
		if (j != -1) {
			base.add(Math.min(indice, j));
			base.add(Math.max(indice, j));
		}

		return base;

	}

	public ArrayList<Integer> findLoopBackward(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		if (indice < 0 || indice >= get_listeBases().size())
			return base;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if (j != -1) {
			i = Math.min(i, j) - 1;
			indice = i;
		}
		if (i < 0) {
			return base;
		}
		do {
			base.add(i);
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
			} else {
				i = (i + get_listeBases().size() - 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return base;
	}

	public ArrayList<Integer> findHelix(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		if (get_listeBases().get(indice).getElementStructure() != -1) {
			list.add(indice);
			list.add(get_listeBases().get(indice).getElementStructure());
			int i = 1, prec = get_listeBases().get(indice)
					.getElementStructure();
			while (indice + i < get_listeBases().size()
					&& get_listeBases().get(indice + i).getElementStructure() != -1
					&& get_listeBases().get(indice + i).getElementStructure() == prec - 1) {
				list.add(indice + i);
				list.add(get_listeBases().get(indice + i)
						.getElementStructure());
				prec = get_listeBases().get(indice + i).getElementStructure();
				i++;
			}
			i = -1;
			prec = get_listeBases().get(indice).getElementStructure();
			while (indice + i >= 0
					&& get_listeBases().get(indice + i).getElementStructure() != -1
					&& get_listeBases().get(indice + i).getElementStructure() == prec + 1) {
				list.add(indice + i);
				list.add(get_listeBases().get(indice + i)
						.getElementStructure());
				prec = get_listeBases().get(indice + i).getElementStructure();
				i--;
			}
		}
		return list;
	}

	public ArrayList<Integer> find3Prime(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		boolean over = false;
		while ((indice >= 0) && !over) {
			over = (get_listeBases().get(indice).getElementStructure() != -1);
			indice--;
		}
		indice++;
		if (over) {
			indice++;
		}
		for (int i = indice; i < get_listeBases().size(); i++) {
			list.add(i);
			if (get_listeBases().get(i).getElementStructure() != -1) {
				return new ArrayList<Integer>();
			}
		}
		return list;
	}

	public ArrayList<Integer> find5Prime(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i <= indice; i++) {
			list.add(i);
			if (get_listeBases().get(i).getElementStructure() != -1) {
				return new ArrayList<Integer>();
			}
		}
		return list;
	}

	public Double get_spaceBetweenBases() {
		return _spaceBetweenBases;
	}

	public void set_spaceBetweenBases(Double betweenBases) {
		_spaceBetweenBases = betweenBases;
	}

	public static Double angle(Point2D.Double p1, Point2D.Double p2,
			Point2D.Double p3) {
		Double alpha = Math.atan2(p1.y - p2.y, p1.x - p2.x);
		Double beta = Math.atan2(p3.y - p2.y, p3.x - p2.x);
		Double angle = (beta - alpha);

		// Correction de l'angle pour le resituer entre 0 et 2PI
		while (angle < 0.0 || angle > 2 * Math.PI) {
			if (angle < 0.0)
				angle += 2 * Math.PI;
			else if (angle > 2 * Math.PI)
				angle -= 2 * Math.PI;
		}
		return angle;
	}

	public ArrayList<Integer> findNonPairedBaseGroup(Integer get_nearestBase) {
		// detection 3', 5', bulge
		ArrayList<Integer> list = new ArrayList<Integer>();
		int indice = get_nearestBase;
		boolean nonpairedUp = true, nonpairedDown = true;
		while (indice < get_listeBases().size() && nonpairedUp) {
			if (get_listeBases().get(indice).getElementStructure() == -1) {
				list.add(indice);
				indice++;
			} else {
				nonpairedUp = false;
			}
		}
		indice = get_nearestBase - 1;
		while (indice >= 0 && nonpairedDown) {
			if (get_listeBases().get(indice).getElementStructure() == -1) {
				list.add(indice);
				indice--;
			} else {
				nonpairedDown = false;
			}
		}
		return list;
	}

	public boolean getDrawn() {
		return _drawn;
	}

	public ArrayList<ModeleBP> getStructureAux() {
		return _structureAux;
	}

	public int getIndexFromBaseNumber(int num) {
		for (int i = 0; i < this._listeBases.size(); i++) {
			if (_listeBases.get(i).getBaseNumber() == num) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Adds a base pair to this RNA's structure. Tries to add it to the
	 * secondary structure first, eventually adding it to the 'tertiary'
	 * interactions if it clashes with the current secondary structure.
	 * 
	 * @param baseNumber5
	 *            - Base number of the origin of this base pair
	 * @param baseNumber3
	 *            - Base number of the destination of this base pair
	 */

	public void addBPToStructureUsingNumbers(int baseNumber5, int baseNumber3) {
		int i = getIndexFromBaseNumber(baseNumber5);
		int j = getIndexFromBaseNumber(baseNumber3);
		addBPToStructure(i, j);
	}

	/**
	 * Adds a base pair to this RNA's structure. Tries to add it to the
	 * secondary structure first, possibly adding it to the 'tertiary'
	 * interactions if it clashes with the current secondary structure.
	 * 
	 * @param number5
	 *            - Base number of the origin of this base pair
	 * @param number3
	 *            - Base number of the destination of this base pair
	 */

	public void addBPToStructureUsingNumbers(int number5, int number3, ModeleBP msbp) {
		addBPToStructure(getIndexFromBaseNumber(number5), getIndexFromBaseNumber(number3), msbp);
	}
	
	public void addBPToStructure(int index5, int index3) {
		int i = index5;
		int j = index3;
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleBP msbp = new ModeleBP(part5, part3);
		addBPToStructure(i, j, msbp);
	}

	public void addBPToStructure(int index5, int index3, ModeleBP msbp)
	{
		int i = index5;
		int j = index3;
		
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		if (i != -1) {
			if ((_listeBases.get(i).getElementStructure() != -1)
					|| (_listeBases.get(j).getElementStructure() != -1)) {
				addBPAux(i, j, msbp);
				return;
			}
			for (int k = i + 1; k < j; k++) {
				ModeleBase tmp = _listeBases.get(k);
				int l = tmp.getElementStructure();
				if (l != -1) {
					if ((l <= i) || (l >= j)) {
						addBPAux(i, j, msbp);
						return;
					}
				}
			}
			addBP(i, j, msbp);
		}
	}
	
	public void removeBP(ModeleBP ms)
	{
 		if (_structureAux.contains(ms))
 		{ _structureAux.remove(ms); }
 		else
 		{
 			ModeleBase m5 = ms.getPartner5();
 			ModeleBase m3 = ms.getPartner3();
 			int i = m5.getIndex();
 			int j = m3.getIndex();
 			if ((m5.getElementStructure()==m3.getIndex())&&(m3.getElementStructure()==m5.getIndex()))
 			{
 				m5.removeElementStructure();
 				m3.removeElementStructure();
 			}
 		}
	}
	

	public void addBP(int i, int j) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}

		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleBP msbp = new ModeleBP(part5, part3);
		addBP(i, j, msbp);
	}

	public void addBP(int i, int j, ModeleBP msbp) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		msbp.setPartner5(part5);
		msbp.setPartner3(part3);
		part5.setElementStructure(j, msbp);
		part3.setElementStructure(i, msbp);
	}

	public void addBPAux(int i, int j) {
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleBP msbp = new ModeleBP(part5, part3);
		addBPAux(i, j, msbp);
	}

	public void addBPAux(int i, int j, ModeleBP msbp) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		msbp.setPartner5(part5);
		msbp.setPartner3(part3);
		_structureAux.add(msbp);
	}


	public ArrayList<ModeleBP> getBPsAt(int i) {
		ArrayList<ModeleBP> result = new ArrayList<ModeleBP>();
		if (_listeBases.get(i).getElementStructure() != -1) {
			result.add(_listeBases.get(i).getStyleBP());
		}
		for (int k = 0; k < _structureAux.size(); k++) {
			ModeleBP bp = _structureAux.get(k);
			if ((bp.getPartner5().getIndex() == i)
					|| (bp.getPartner3().getIndex() == i)) {
				result.add(bp);
			}
		}
		return result;
		
	}

	public ModeleBP getBPStyle(int i, int j) {
		ModeleBP result = null;
		if (i > j) {
			int k = j;
			j = i;
			i = k;
		}
		if (_listeBases.get(i).getElementStructure() == j) {
			result = _listeBases.get(i).getStyleBP();
		}
		for (int k = 0; k < _structureAux.size(); k++) {
			ModeleBP bp = _structureAux.get(k);
			if ((bp.getPartner5().getIndex() == i)
					&& (bp.getPartner3().getIndex() == j)) {
				result = bp;
			}
		}
		return result;
	}

	public ArrayList<ModeleBP> getSecStrBPs()
	{
		ArrayList<ModeleBP> result = new ArrayList<ModeleBP>(); 
		for (int i=0;i<this.getSize();i++) 
		{
			ModeleBase mb = _listeBases.get(i); 
			int k = mb.getElementStructure();
			if ((k!=-1) && (k>i))
			{  result.add(mb.getStyleBP()) ;  }
		}
		return result;
	}

	public ArrayList<ModeleBP> getAuxBPs()
	{
		ArrayList<ModeleBP> result = new ArrayList<ModeleBP>(); 
		for (ModeleBP bp : _structureAux) 
		{
			result.add(bp);
		}
		return result;
	}

	public ArrayList<ModeleBP> getAllBPs()
	{
		ArrayList<ModeleBP> result = new ArrayList<ModeleBP>();
		result.addAll(getSecStrBPs());
		result.addAll(getAuxBPs());
		return result;
	}

	
	public ArrayList<ModeleBP> getAuxBPs(int i)
	{
		ArrayList<ModeleBP> result = new ArrayList<ModeleBP>(); 
		for (ModeleBP bp : _structureAux) 
		{
			if ((bp.getPartner5().getIndex() == i) || (bp.getPartner3().getIndex() == i)) {
				result.add(bp);
			}
		}
		return result;
	}
	
	public void setBaseInnerColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_inner_color(c);
		}
	}

	public void setBaseNumbersColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_number_color(c);
		}
	}

	public void setBaseNameColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_name_color(c);
		}
	}

	public void setBaseOutlineColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_outline_color(c);
		}
	}

	public String getName()
	{
		return _name;
	}

	public void setName(String n)
	{
		_name = n;
	}
	
	public ArrayList<TextAnnotation> getAnnotations()
	{
		return _listeAnnotations;
	}

	public boolean removeAnnotation(TextAnnotation t)
	{
		return _listeAnnotations.remove(t);
	}

	public void addAnnotation(TextAnnotation t)
	{
		_listeAnnotations.add(t);
	}

	public void removeAnnotation(String filter)
	{
		ArrayList<TextAnnotation> condamne = new ArrayList<TextAnnotation>();
		for (TextAnnotation t : _listeAnnotations)
		{
			if (t.getTexte().contains(filter))
			{
				condamne.add(t);				
			}
		}
		for (TextAnnotation t : condamne)
		{
			_listeAnnotations.remove(t);
		}
	}

	public void clearAnnotations()
	{
		_listeAnnotations.clear();
	}
	
	private boolean _strandEndsAnnotated = false;
	
	public void autoAnnotateStrandEnds()
	{
		if (! _strandEndsAnnotated)
		{
		int tailleListBases=_listeBases.size();
		boolean endAnnotate =false;
		addAnnotation(new TextAnnotation("5'", _listeBases.get(0)));
		for(int i=0; i<_listeBases.size()-1; i++){
			int realposA = _listeBases.get(i).getBaseNumber();
			int realposB = _listeBases.get(i+1).getBaseNumber();
			if(realposB-realposA!=1){
				addAnnotation(new TextAnnotation("3'", _listeBases.get(i)));
				addAnnotation(new TextAnnotation("5'", _listeBases.get(i+1)));
				if(i+1==_listeBases.size()-1){
					endAnnotate=true;
				}
			}
		}
		if(!endAnnotate){
			addAnnotation(new TextAnnotation("3'", _listeBases.get(tailleListBases-1)));
		}
		_strandEndsAnnotated = true;
		}
		else
		{
			removeAnnotation("3'");
			removeAnnotation("5'");
			_strandEndsAnnotated = false;
		}
	}
	
	public void autoAnnotateHelices()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbH = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
				ModeleBase mb = _listeBases.get(i); 
				int j = mb.getElementStructure(); 
				if (j==-1)
				{ p.push(i+1);	}
				else
				{
					if (j>i)
					{
						ModeleBase mbp = _listeBases.get(j); 
						p.push(j+1);
						ArrayList<ModeleBase> h = new ArrayList<ModeleBase>(); 
						int k = 1;
						while(mb.getElementStructure()==mbp.getIndex())
						{
							h.add(mb);
							h.add(mbp);
							mb = _listeBases.get(i+k);
							mbp = _listeBases.get(j-k); 
						
							k++;
						}
						try {
							addAnnotation(new TextAnnotation("H"+nbH++, h , TextAnnotation.HELIX));
						} catch (Exception e) {
							e.printStackTrace();
						}
						p.push(i+k);
					}
				}
			}
		}
	}

	public void autoAnnotateTerminalLoops()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbT = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
			ModeleBase mb = _listeBases.get(i); 
			int j = mb.getElementStructure(); 
			if (j==-1)
			{ 	
				int k = 1;
				ArrayList<ModeleBase> t = new ArrayList<ModeleBase>();
				while ((i+k<getSize())&&(mb.getElementStructure()==-1))
				{
					t.add(mb);
					mb = _listeBases.get(i+k);
					k++;
				}
				if (mb.getElementStructure()!=-1)
				{
					if (mb.getElementStructure()==i-1)
					{
						try {
							t.add(_listeBases.get(i-1));
							t.add(_listeBases.get(i+k-1));
							addAnnotation(new TextAnnotation("T"+nbT++, t , TextAnnotation.LOOP));
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
					p.push(i+k-1);
				}
				
			}
			else
			{
				if (j>i)
				{
					p.push(j+1);
					p.push(i+1);
				}
			}
			}
		}
	}

	public void autoAnnotateInteriorLoops()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbT = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
				ModeleBase mb = _listeBases.get(i); 
				int j = mb.getElementStructure(); 
				if (j==-1)
				{ 	
					int k = i+1;
					ArrayList<ModeleBase> t = new ArrayList<ModeleBase>();
					boolean terminal = true;
					while ((k<getSize())&&((mb.getElementStructure()>=i)||(mb.getElementStructure()==-1)))
					{
						t.add(mb);
						mb = _listeBases.get(k);
						if ((mb.getElementStructure()==-1)||(mb.getElementStructure()<k))
							k++;
						else 
						{
							p.push(k);
							terminal = false;
							k = mb.getElementStructure();
						}
					}
					if (mb.getElementStructure()!=-1)
					{
						if ((mb.getElementStructure()==i-1)&& !terminal)
						{
							try {
								t.add(_listeBases.get(i-1));
								t.add(_listeBases.get(k-1));
								addAnnotation(new TextAnnotation("I"+nbT++, t , TextAnnotation.LOOP));
							} catch (Exception e) {
								e.printStackTrace();
							}
							p.push(k-1);
						}
					}	
				}
				else
				{
				if (j>i)
				{
					p.push(i+1);
				}
			}
		}
		}
	}
	
	@SuppressWarnings("unchecked")
	public TextAnnotation getAnnotation(int type, ModeleBase base)
	{
		TextAnnotation result = null;
		for(TextAnnotation t : _listeAnnotations)
		{
			if (t.getType()==type)
			{
				switch(type)
				{
				    case(TextAnnotation.BASE):
				    	if (base == (ModeleBase) t.getAncrage())
				    	  return t;
					break;
				    case(TextAnnotation.HELIX):
				    case(TextAnnotation.LOOP):
				    {
				    	ArrayList<ModeleBase> mbl = (ArrayList<ModeleBase>) t.getAncrage();
				    	if (mbl.contains(base))
				    	  return t;
				    }
					break;
				}
			}
		}
		return result;
	}
	
	private ArrayList<ChemProbAnnotation> _ChemProbAnnotations = new ArrayList<ChemProbAnnotation>();
	
	public void addChemProbAnnotation(ChemProbAnnotation cpa)
	{
		_ChemProbAnnotations.add(cpa);
	}
	
	public ArrayList<ChemProbAnnotation> getChemProbAnnotations()
	{
		return _ChemProbAnnotations;
	}

	public void setColorMapValues(Double[] values, ModeleColorMap cm)
	{
		setColorMapValues(values, cm, false);
	}
	
	public void adaptColorMapToValues(ModeleColorMap cm)
	{
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		for (int i=0;i<Math.min(_listeBases.size(),_listeBases.size());i++)
		{
			ModeleBase mb = _listeBases.get(i);
			max = Math.max(max,mb.getValue());
			min = Math.min(min,mb.getValue());
		}
		cm.rescale(min, max);		
	}
	
	  public void readValues(Reader r, ModeleColorMap cm)
	  {
		  try {
			  StreamTokenizer st = new StreamTokenizer(r);
			  st.eolIsSignificant(true);
			  ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();
			  ArrayList<Double> curVals = new ArrayList<Double>(); 
			  int type = st.nextToken();
			  while (type != StreamTokenizer.TT_EOF)
			  {
				  switch(type)
				  {
				  case (StreamTokenizer.TT_NUMBER): 
					  curVals.add(st.nval);
				  break;
				  case (StreamTokenizer.TT_EOL):
					  if (curVals.size()>0)
					  {
						  vals.add(curVals);
						  curVals = new ArrayList<Double>();
					  }
				  break;
				  }
				  type = st.nextToken();
			  }
			  if (curVals.size()>0) 
				  vals.add(curVals);
			  
			  Double[] v = new Double[vals.size()];
			  for (int i=0;i<Math.min(vals.size(),getSize());i++)
			  {
				  ArrayList<Double> tab = vals.get(i);
				  v[i] = tab.get(tab.size()-1);
			  }
			  setColorMapValues(v, cm,  true);
		  } catch (IOException e) {
			  e.printStackTrace();
		  }
	  }


	public void setColorMapValues(Double[] values, ModeleColorMap cm, boolean rescaleColorMap)
	{
		if (values.length>0)
		{
			for (int i=0;i<Math.min(values.length,_listeBases.size());i++)
			{
				ModeleBase mb = _listeBases.get(i);
				mb.setValue(values[i]);
			}
			if (rescaleColorMap){
				adaptColorMapToValues(cm);
			}
		}
	}

	public Double[] getColorMapValues()
	{
		Double[] values = new Double[_listeBases.size()];
		for (int i=0;i<_listeBases.size();i++)
		{
			values[i] = _listeBases.get(i).getValue();
		}
		return values;
	}
	
	public void rescaleColorMap(ModeleColorMap cm)
	{
		Double max = Double.MIN_VALUE;
		Double min = Double.MAX_VALUE;
		for (int i=0;i<_listeBases.size();i++)
		{
			Double value = _listeBases.get(i).getValue();
			max = Math.max(max, value);
			min = Math.min(min, value);
		}
		cm.rescale(min, max);
	}
	

	public void setSequence(String s)
	{
		setSequence(RNA.explodeSequence(s));
	}
	public void setSequence(List<String> s)
	{
		int i = 0;
		int j = 0;
		while ((i<s.size())&&(j<_listeBases.size()))
		{
			ModeleBase mb = _listeBases.get(j);
			if (mb instanceof ModeleBaseNucleotide)
			{
				((ModeleBaseNucleotide)mb).set_c(s.get(i));
				i++; j++;
			}
			else if (mb instanceof ModeleBasesComparison)
			{
				((ModeleBasesComparison)mb).set_base1(((s.get(i).length()>0)?s.get(i).charAt(0):' '));
				((ModeleBasesComparison)mb).set_base2(((s.get(i+1).length()>0)?s.get(i+1).charAt(0):' '));
				i+=2;j++;
			}
			else j++;
 		}
	}
	
	public void eraseSequence()
	{
		int j = 0;
		while ((j<_listeBases.size()))
		{
			ModeleBase mb = _listeBases.get(j);
			if (mb instanceof ModeleBaseNucleotide)
			{
				((ModeleBaseNucleotide)mb).set_c("");
				j++;
			}
			else if (mb instanceof ModeleBasesComparison)
			{
				((ModeleBasesComparison)mb).set_base1(' ');
				((ModeleBasesComparison)mb).set_base2(' ');
				j++;
			}
			else j++;
 		}		
	}
	
	
    public RNA clone ()
    {
        try
        {
            ByteArrayOutputStream out = new ByteArrayOutputStream ();
            ObjectOutputStream oout = new ObjectOutputStream (out);
            oout.writeObject (this);
            
            ObjectInputStream in = new ObjectInputStream (
                new ByteArrayInputStream (out.toByteArray ()));
            return (RNA)in.readObject ();
        }
        catch (Exception e)
        {
            throw new RuntimeException ("cannot clone class [" +
                this.getClass ().getName () + "] via serialization: " +
                e.toString ());
        }
    }
    
    public ModeleBase getBaseAt(int index)
    {
    	return this._listeBases.get(index);
    }

    public ArrayList<ModeleBase> getBasesAt(Collection<? extends Integer> indices)
    {
    	ArrayList<ModeleBase> mbs = new ArrayList<ModeleBase>();
    	for (int i: indices)
    	{
    		mbs.add(getBaseAt(i));
    	}
    	return mbs;
    }

    public ArrayList<ModeleBase> getBasesBetween(int from, int to)
    {
    	ArrayList<ModeleBase> mbs = new ArrayList<ModeleBase>();
    	int bck = Math.min(from, to);
    	to = Math.max(from, to);
    	from = bck;
    	for (int i=from;i<=to;i++)
    	{  mbs.add(getBaseAt(i)); }
    	return mbs;
    }

    public void addHighlightRegion(HighlightRegionAnnotation n)
    {
    	_listeRegionHighlights.add(n);    	
    }

	public void removeHighlightRegion(HighlightRegionAnnotation n)
	{
		_listeRegionHighlights.remove(n);
	}
   
	public void removeChemProbAnnotation(ChemProbAnnotation a)
	{
		_ChemProbAnnotations.remove(a);
	}
	
    public void addHighlightRegion(int from, int to, Color fill, Color outline, double radius)
    {
    	_listeRegionHighlights.add(new HighlightRegionAnnotation(getBasesBetween(from, to),fill,outline,radius));    	
    }
    
    
    public void addHighlightRegion(int from, int to)
    {
    	_listeRegionHighlights.add(new HighlightRegionAnnotation(getBasesBetween(from, to)));    	
    }

    public ArrayList<HighlightRegionAnnotation> getHighlightRegion()
    {
    	return _listeRegionHighlights;    	
    }


	/**
	 * Rotates the RNA coordinates by a certain angle
	 * 
	 * @param angleDegres
	 *            Rotation angle, in degrees
	 */
	public void globalRotation(Double angleDegres) {
		if (_listeBases.size() > 0) {

			// angle en radian
			Double angle = angleDegres * Math.PI / 180;

			// initialisation du minimum et dumaximum
			Double maxX = _listeBases.get(0).getCoords().x;
			Double maxY = _listeBases.get(0).getCoords().y;
			Double minX = _listeBases.get(0).getCoords().x;
			Double minY = _listeBases.get(0).getCoords().y;
			// mise a jour du minimum et du maximum
			for (int i = 0; i < _listeBases.size(); i++) {
				if (_listeBases.get(i).getCoords().getX() < minX)
					minX = _listeBases.get(i).getCoords().getX();
				if (_listeBases.get(i).getCoords().getY() < minY)
					minY = _listeBases.get(i).getCoords().getY();
				if (_listeBases.get(i).getCoords().getX() > maxX)
					maxX = _listeBases.get(i).getCoords().getX();
				if (_listeBases.get(i).getCoords().getX() > maxY)
					maxY = _listeBases.get(i).getCoords().getY();
			}
			// creation du point central
			Point2D.Double centre = new Point2D.Double((maxX - minX) / 2,
					(maxY - minY) / 2);
			Double x, y;
			for (int i = 0; i < _listeBases.size(); i++) {
				// application de la rotation au centre de chaque base
				// x' = cos(theta)*(x-xc) - sin(theta)*(y-yc) + xc
				x = Math.cos(angle)
						* (_listeBases.get(i).getCenter().getX() - centre.x)
						- Math.sin(angle)
						* (_listeBases.get(i).getCenter().getY() - centre.y)
						+ centre.x;
				// y' = sin(theta)*(x-xc) + cos(theta)*(y-yc) + yc
				y = Math.sin(angle)
						* (_listeBases.get(i).getCenter().getX() - centre.x)
						+ Math.cos(angle)
						* (_listeBases.get(i).getCenter().getY() - centre.y)
						+ centre.y;
				_listeBases.get(i).setCenter(
						new Point2D.Double(x, y));

				// application de la rotation au coordonnees de chaque
				// base
				// x' = cos(theta)*(x-xc) - sin(theta)*(y-yc) + xc
				x = Math.cos(angle)
						* (_listeBases.get(i).getCoords().getX() - centre.x)
						- Math.sin(angle)
						* (_listeBases.get(i).getCoords().getY() - centre.y)
						+ centre.x;
				// y' = sin(theta)*(x-xc) + cos(theta)*(y-yc) + yc
				y = Math.sin(angle)
						* (_listeBases.get(i).getCoords().getX() - centre.x)
						+ Math.cos(angle)
						* (_listeBases.get(i).getCoords().getY() - centre.y)
						+ centre.y;
				_listeBases.get(i).setCoords(
						new Point2D.Double(x, y));
			}
		}
	}
	public boolean testDirectionality(int i, int j, int k) {

		// Which direction are we heading toward?
		Point2D.Double pi = getCoords(i);
		Point2D.Double pj = getCoords(j);
		Point2D.Double pk = getCoords(k);
		return testDirectionality(pi, pj, pk);
	}

	public static boolean testDirectionality(Point2D.Double pi, Point2D.Double pj, Point2D.Double pk) {

		// Which direction are we heading toward?
		double test = (pj.x - pi.x) * (pk.y - pj.y) - (pj.y - pi.y)
				* (pk.x - pj.x);
		return test < 0.0;
	}
	public double getOrientation()
	{
		double maxDist = Double.MIN_VALUE;
		double angle = 0;
		for(int i=0; i< _listeBases.size();i++)
		{
			ModeleBase b1 = _listeBases.get(i);
			for(int j=i+1; j< _listeBases.size();j++)
			{
				ModeleBase b2 = _listeBases.get(j);
				Point2D.Double p1 = b1._coords.toPoint2D();
				Point2D.Double p2 = b2._coords.toPoint2D();
				double dist = p1.distance(p2); 
				if(dist>maxDist)
				{
					maxDist = dist;
					angle = VueUI.computeAngle(p1, p2);
				}
			}
		}
		return angle;
	}
	
	public boolean hasVirtualLoops()
	{
		boolean consecutiveBPs = false;
		for(int i=0;i<_listeBases.size();i++)
		{
			int j =_listeBases.get(i).getElementStructure(); 
			if (j==i+1)
			{ consecutiveBPs = true; }
			
		}
		return ((_drawMode!=DRAW_MODE_LINEAR)
			  &&(_drawMode!=DRAW_MODE_CIRCULAR)
			  &&(consecutiveBPs));
	}
}


