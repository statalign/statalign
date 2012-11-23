package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;

import javax.swing.undo.UndoManager;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.controlers.ControleurBlinkingThread;
import fr.orsay.lri.varna.controlers.ControleurClicMovement;
import fr.orsay.lri.varna.controlers.ControleurDraggedMolette;
import fr.orsay.lri.varna.controlers.ControleurInterpolator;
import fr.orsay.lri.varna.controlers.ControleurMolette;
import fr.orsay.lri.varna.controlers.ControleurVARNAPanelKeys;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.export.VueVARNAGraphics;
import fr.orsay.lri.varna.models.rna.ModeleBP;
import fr.orsay.lri.varna.models.rna.RNA;
import fr.orsay.lri.varna.views.VueUI;

public class StructureGUI extends VARNAPanel{

	/**
	 * This class implements a graphical interface for displaying the current consensus structure. RNA
	 * graphics are provided by VARNA.
	 * 
	 * @author Preeti Arunapuram
	 */
	private static final long serialVersionUID = 1L;
	private static String sequence;
	private static String structure;
	
	private Graphics2D g2;
	private ColorGradient cg = new ColorGradient(Color.WHITE, Color.BLACK);
	
	private RNA _RNA = new RNA();
	
	
	public String title;
	public boolean probMode = true;
	
	private double _scaleFactor = 1.0;
	
	private float[][] probs;
	

	private VARNAConfig _conf = new VARNAConfig();
	
	UndoManager _manager;
	
	private ControleurBlinkingThread _blink;

	private Point _translation;
	private boolean _horsCadre;
	private boolean _premierAffichage;


	private ControleurInterpolator _interpolator;

	private VueUI _UI = new VueUI(this);

	
	public StructureGUI(String title) throws ExceptionNonEqualLength {
		sequence = null;
		structure = null;
		this.title = title;
		
		
	}
	
	public void paintComponent(Graphics g) {
		//setBackground(Color.PINK);
		//setBaseInnerColor(Color.GREEN);
		//setBaseNameColor(Color.BLUE);
		
		super.paintComponent(g);
		
		g2 = (Graphics2D)g;
		
		if(structure == null || sequence == null) {
			g2.drawString("Waiting for data..", 30, 30);
            return;
		}
	
		g2.setPaint(Color.BLACK);
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		g2.drawString(title, 10, 20);
		
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 12));
		
		if(!probMode) g2.drawString("Normal Mode", 10, 40);
		else g2.drawString("Probability Mode", 10, 40);
		
		g2.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		
	}
	
	public void updateAndDraw(String seq, String str) {		
		sequence = seq;
		structure = str;
		
		repaint();
		
		super.drawRNAInterpolated(sequence, structure, RNA.DRAW_MODE_RADIATE);
		setBackground(VARNAConfig.DEFAULT_BACKGROUND_COLOR);
		_manager = new UndoManager();
		_manager.setLimit(10000);
		_UI.addUndoableEditListener(_manager);
		

		_blink = new ControleurBlinkingThread(this,
				ControleurBlinkingThread.DEFAULT_FREQUENCY, 0, 1.0, 0.0, 0.2);
		_blink.start();

		_premierAffichage = true;
		_translation = new Point(0, 0);

		_horsCadre = false;
		this.setFont(_conf._fontBasesGeneral);

		// ajout des controleurs au VARNAPanel
		ControleurClicMovement controleurClicMovement = new ControleurClicMovement(	this);
		this.addMouseListener(controleurClicMovement);
		this.addMouseMotionListener(controleurClicMovement);
		this.addMouseWheelListener(new ControleurMolette(this));

		ControleurDraggedMolette ctrlDraggedMolette = new ControleurDraggedMolette(
				this);
		this.addMouseMotionListener(ctrlDraggedMolette);
		this.addMouseListener(ctrlDraggedMolette);

		ControleurVARNAPanelKeys ctrlKey = new ControleurVARNAPanelKeys(this);
		this.addKeyListener(ctrlKey);
		this.addFocusListener(ctrlKey);
		
		_interpolator = new ControleurInterpolator(this);
		_interpolator.start();
	}
	
	//@Override
	public void drawSymbol(VueVARNAGraphics g2D,double posx, double posy,
			double normx, double normy, double radius, boolean isCIS,
			ModeleBP.Edge e) {

		switch (e) {
		case WATSON_CRICK:
			if (isCIS) {
				g2D.fillCircle( (posx - (radius) / 2.0),
						(posy - (radius) / 2.0),  radius);
			} else {
				Color bck = g2D.getColor();
				g2D.setColor(Color.white);
				g2D.fillCircle( posx - (radius) / 2.0,
						 (posy - (radius) / 2.0),  (radius));
				g2D.setColor(bck);
				g2D.drawCircle( (posx - (radius) / 2.0),
						 (posy - (radius) / 2.0),  (radius));
			}
			break;
		case HOOGSTEEN: {
			GeneralPath p2 = new GeneralPath();
			p2.moveTo((float) (posx - radius * normx / 2.0 - radius * normy
					/ 2.0), (float) (posy - radius * normy / 2.0 + radius
					* normx / 2.0));
			p2.lineTo((float) (posx + radius * normx / 2.0 - radius * normy
					/ 2.0), (float) (posy + radius * normy / 2.0 + radius
					* normx / 2.0));
			p2.lineTo((float) (posx + radius * normx / 2.0 + radius * normy
					/ 2.0), (float) (posy + radius * normy / 2.0 - radius
					* normx / 2.0));
			p2.lineTo((float) (posx - radius * normx / 2.0 + radius * normy
					/ 2.0), (float) (posy - radius * normy / 2.0 - radius
					* normx / 2.0));
			p2.closePath();

			if (isCIS) {
				g2D.fill(p2);
			} else {
				Color bck = g2D.getColor();
				g2D.setColor(Color.white);
				g2D.fill(p2);
				g2D.setColor(bck);
				g2D.draw(p2);
			}
		}
			break;
		case SUGAR: {
			double ix = radius * normx / 2.0;
			double iy = radius * normy / 2.0;
			double jx = radius * normy / 2.0;
			double jy = -radius * normx / 2.0;

			GeneralPath p2 = new GeneralPath();
			p2.moveTo((float) (posx - ix + jx), (float) (posy - iy + jy));
			p2.lineTo((float) (posx + ix + jx), (float) (posy + iy + jy));
			p2.lineTo((float) (posx - jx), (float) (posy - jy));
			p2.closePath();

			if (isCIS) {
				g2D.fill(p2);
			} else {
				Color bck = g2D.getColor();
				g2D.setColor(Color.white);
				g2D.fill(p2);
				g2D.setColor(bck);
				g2D.draw(p2);
			}
		}
			break;
		}
	}
	
	//@Override
	public void drawBasePair(VueVARNAGraphics g2D,Point2D.Double orig,
			Point2D.Double dest, ModeleBP style, double newRadius) {
		
		double dx = dest.x - orig.x;
		double dy = dest.y - orig.y;
		double dist = Math.sqrt((dest.x - orig.x) * (dest.x - orig.x)
				+ (dest.y - orig.y) * (dest.y - orig.y));
		dx /= dist;
		dy /= dist;
		double nx = -dy;
		double ny = dx;
		orig = new Point2D.Double(orig.x+newRadius*dx,orig.y+newRadius*dy);
		dest = new Point2D.Double(dest.x-newRadius*dx,dest.y-newRadius*dy);
		if (_conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_LW) 
		{
			double radiusCircle = ((RNA.BASE_PAIR_DISTANCE - _RNA.BASE_RADIUS) / 5.0)
					* this._scaleFactor;
			if (style.isCanonical()) 
			{
				g2D.setStrokeThickness(probs[style.getIndex3()][style.getIndex5()]);
				((Graphics2D) g2D).setPaint(cg.getColor(probs[style.getIndex3()][style.getIndex5()]));
				if (style.isCanonicalGC()) 
				{
					if ((orig.x != dest.x) || (orig.y != dest.y)) {
						nx *= this._scaleFactor * _RNA.BASE_RADIUS / 4.0;
						ny *= this._scaleFactor * _RNA.BASE_RADIUS / 4.0;
						g2D.drawLine( (orig.x + nx), (orig.y + ny),
								 (dest.x + nx), (dest.y + ny));
						g2D.drawLine( (orig.x - nx),  (orig.y - ny),
								 (dest.x - nx),  (dest.y - ny));
					}
				} 
				else if (style.isCanonicalAU()){
					g2D.drawLine(orig.x, orig.y, dest.x,
							 dest.y);
				}
				else if (style.isWobbleUG()){
					double cx = (dest.x + orig.x) / 2.0;
					double cy = (dest.y + orig.y) / 2.0;
					g2D.drawLine( orig.x,  orig.y,  dest.x,
							dest.y);
					drawSymbol(g2D, cx, cy, nx, ny, radiusCircle, false, ModeleBP.Edge.WATSON_CRICK);
				}
				else  
				{
					double cx = (dest.x + orig.x) / 2.0;
					double cy = (dest.y + orig.y) / 2.0;
					g2D.drawLine( orig.x,  orig.y,  dest.x,
							dest.y);
					drawSymbol(g2D, cx, cy, nx, ny, radiusCircle, style
							.isCIS(), style.getEdgePartner5());
				} 
			} else {
				ModeleBP.Edge p1 = style.getEdgePartner5();
				ModeleBP.Edge p2 = style.getEdgePartner3();
				double cx = (dest.x + orig.x) / 2.0;
				double cy = (dest.y + orig.y) / 2.0;
				g2D.drawLine( orig.x, orig.y,dest.x,
						dest.y);
				if (p1 == p2) {
					drawSymbol(g2D, cx, cy, nx, ny, radiusCircle, style
							.isCIS(), p1);

				} else {
					double vdx = (dest.x - orig.x);
					double vdy = (dest.y - orig.y);
					vdx /= 6.0;
					vdy /= 6.0;
					drawSymbol(g2D, cx + vdx, cy + vdy, nx, ny, radiusCircle,
							style.isCIS(), p2);
					drawSymbol(g2D, cx - vdx, cy - vdy, nx, ny, radiusCircle,
							style.isCIS(), p1);
				}
			}
		} else if (_conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_SIMPLE) {
			g2D.drawLine( orig.x,  orig.y,dest.x, dest.y);
		} else if (_conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_RNAVIZ) {
			double xcenter = (orig.x + dest.x) / 2.0;
			double ycenter = (orig.y + dest.y) / 2.0;
			double radius = Math.max(4.0 * this._scaleFactor, 1.0);
			g2D.fillCircle( (xcenter - radius), (ycenter - radius),  (2.0 * radius));
		}
	}
	
	public void setMatrix(float[][] probMatrix) {
		probs = probMatrix;
	}
	
}
	
	
	
	
	
	
	

