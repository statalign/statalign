package fr.orsay.lri.varna.controlers;

import java.awt.Color;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.models.rna.ModeleColorMap;
import fr.orsay.lri.varna.models.rna.RNA;

public class ControleurScriptParser {
    private static String SCRIPT_ERROR_PREFIX = "Error";	


	private static class Command{
		Function _f;
		Vector<Argument> _argv; 
		public Command(Function f, Vector<Argument> argv)
		{
			_f = f;
			_argv = argv;
		}
	}; 

	private static abstract class Argument{
		ArgumentType _t;
		
		public Argument(ArgumentType t)
		{
			_t = t;
		}
		
		public ArgumentType getType()
		{
			return _t;
		}
		
		public abstract String toString();
		
	}; 
	

	private static class NumberArgument extends Argument{
		Number _val; 
		public NumberArgument(Number val)
		{
			super(ArgumentType.NUMBER_TYPE);
			_val = val;
		}
		public Number getNumber()
		{
			return _val;
		}
		public String toString()
		{
			return _val.toString();
		}
	}; 
	
	private static class ColorArgument extends Argument{
		Color _val; 
		public ColorArgument(Color val)
		{
			super(ArgumentType.COLOR_TYPE);
			_val = val;
		}
		public Color getColor()
		{
			return _val;
		}
		public String toString()
		{
			return _val.toString();
		}
	}; 

	private static class StringArgument extends Argument{
		String _val; 
		public StringArgument(String val)
		{
			super(ArgumentType.STRING_TYPE);
			_val = val;
		}
		public String toString()
		{
			return _val.toString();
		}
	}; 

	private static class ArrayArgument extends Argument{
		Vector<Argument> _val; 
		public ArrayArgument(Vector<Argument> val)
		{
			super(ArgumentType.ARRAY_TYPE);
			_val = val;
		}
		public int getSize()
		{
			return _val.size();
		}
		public Argument getArgument(int i)
		{
			return _val.get(i);
		}
		public String toString()
		{
			return _val.toString();
		}
	}; 
	
	
	private enum ArgumentType{
		STRING_TYPE,
		NUMBER_TYPE,
		ARRAY_TYPE,
		COLOR_TYPE
	};
    
    
    private enum Function{
    	ERASE_SEQ,
		SET_COLOR_MAP_MIN,
		SET_COLOR_MAP_MAX,
		SET_COLOR_MAP,
		SET_CUSTOM_COLOR_MAP,
		SET_SEQ,
		SET_STRUCT,
		SET_STRUCT_SMOOTH,
		SET_TITLE,
		SET_RNA,
		SET_RNA_SMOOTH,
		SET_VALUES,
		TOGGLE_SHOW_COLOR_MAP,
		REDRAW,
		UNKNOWN,
	};
	
	private static Hashtable<String,Function> _name2Fun = new Hashtable<String,Function>();  
	private static Hashtable<Function,ArgumentType[]> _fun2Prot = new Hashtable<Function,ArgumentType[]>();
	
	
	private static void initFunctions()
	{
		if (_name2Fun.size()>0)
		{ return; }
		{
			String funtxt =        "eraseseq";
			Function fun =         Function.ERASE_SEQ;
			ArgumentType[] proto = {};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "settitle";
			Function fun =         Function.SET_TITLE;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "redraw";
			Function fun =         Function.REDRAW;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setstruct";
			Function fun =         Function.SET_STRUCT;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setseq";
			Function fun =         Function.SET_SEQ;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setrna";
			Function fun =         Function.SET_RNA;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE,ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setstructsmooth";
			Function fun =         Function.SET_STRUCT_SMOOTH;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setrnasmooth";
			Function fun =         Function.SET_RNA_SMOOTH;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE,ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setvalues";
			Function fun =         Function.SET_VALUES;
			ArgumentType[] proto = {ArgumentType.ARRAY_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setcolormap";
			Function fun =         Function.SET_COLOR_MAP;
			ArgumentType[] proto = {ArgumentType.STRING_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setcolormapminvalue";
			Function fun =         Function.SET_COLOR_MAP_MIN;
			ArgumentType[] proto = {ArgumentType.NUMBER_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setcolormapmaxvalue";
			Function fun =         Function.SET_COLOR_MAP_MAX;
			ArgumentType[] proto = {ArgumentType.NUMBER_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "setcustomcolormap";
			Function fun =         Function.SET_CUSTOM_COLOR_MAP;
			ArgumentType[] proto = {ArgumentType.ARRAY_TYPE};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
		{
			String funtxt =        "toggleshowcolormap";
			Function fun =         Function.TOGGLE_SHOW_COLOR_MAP;
			ArgumentType[] proto = {};
			_name2Fun.put(funtxt,fun);_fun2Prot.put(fun,proto);
		}
	}
	
	private static Function getFunction(String f)
	{
		String s = f.trim().toLowerCase();
		if (_name2Fun.containsKey(s))
			return _name2Fun.get(s);
		return Function.UNKNOWN;
	}

	private static ArgumentType[] getPrototype(Function f)
	{
		if (_fun2Prot.containsKey(f))
			return _fun2Prot.get(f);
		return new ArgumentType[0];
	}	
	
    public static void executeScript(VARNAPanel vp, String cmdtxt) throws Exception
    {
    	Vector<Command> cmds = parseScript(cmdtxt);
    	for(int i=0;i<cmds.size();i++)
    	{
    		Command cmd = cmds.get(i);
    		switch(cmd._f)
    		{
    			case ERASE_SEQ:
    			{
    				vp.eraseSequence();
    			}
    			break;
    			case SET_COLOR_MAP_MIN:
    			{
    				vp.setColorMapMinValue(((NumberArgument) cmd._argv.get(0)).getNumber().doubleValue());
    			}
    			break;
    			case SET_COLOR_MAP_MAX:
    			{
    				vp.setColorMapMaxValue(((NumberArgument) cmd._argv.get(0)).getNumber().doubleValue());
    			}
    			break;
    			case SET_COLOR_MAP:
    			{
    				vp.setColorMap(ModeleColorMap.parseColorMap(cmd._argv.get(0).toString()));
    			}
    			break;
    			case SET_CUSTOM_COLOR_MAP:
    			{
    				ModeleColorMap cm = new ModeleColorMap();
    				//System.out.println("a"+cmd._argv.get(0));
    				ArrayArgument arg = (ArrayArgument) cmd._argv.get(0);
    				for (int j=0;j<arg.getSize();j++)
    				{
    					Argument a = arg.getArgument(j);
    					if (a._t==ArgumentType.ARRAY_TYPE)
    					{ 
    	    				//System.out.println("%");
    						ArrayArgument aarg = (ArrayArgument) a; 
    						if (aarg.getSize()==2)
    						{
    							Argument a1 = aarg.getArgument(0);
    							Argument a2 = aarg.getArgument(1);
        	    				//System.out.println("& |"+a1+"| ["+a1.getType()+"] |"+a2+"| ["+a2.getType()+"]");
    							if ((a1.getType()==ArgumentType.NUMBER_TYPE)&&(a2.getType()==ArgumentType.COLOR_TYPE))
    							{
            	    				//System.out.println("+");
    								cm.addColor(((NumberArgument)a1).getNumber().doubleValue(),((ColorArgument)a2).getColor());
    							}
    						}
    					}
    				}    				
    				vp.setColorMap(cm);
    			}
    			break;
    			case SET_TITLE:
    			{
    				vp.setTitle(cmd._argv.get(0).toString());
    			}
    			break;
    			case SET_STRUCT:
    			{
    				String seq = vp.getRNA().getSeq();
    				String str = cmd._argv.get(0).toString();
    				vp.drawRNA(seq, str);
    			}
    			break;
    			case SET_SEQ:
    			{
    				String seq = cmd._argv.get(0).toString();
    				vp.setSequence(seq);
    			}
    			break;
    			case SET_RNA:
    			{
    				String seq = cmd._argv.get(0).toString();
    				String str = cmd._argv.get(1).toString();
    				vp.drawRNA(seq, str);
    			}
    			break;
    			case SET_STRUCT_SMOOTH:
    			{
    				String seq = vp.getRNA().getSeq();
    				String str = cmd._argv.get(0).toString();
    				vp.drawRNAInterpolated(seq, str);
    				vp.repaint();
    			}
    			break;
    			case SET_RNA_SMOOTH:
    			{
    				String seq = cmd._argv.get(0).toString();
    				String str = cmd._argv.get(1).toString();
    				vp.drawRNAInterpolated(seq, str);
    				vp.repaint();
    			}
    			break;
    			case SET_VALUES:
    			{
    				ArrayArgument arg = (ArrayArgument) cmd._argv.get(0);
    				Double[] vals = new Double[arg.getSize()];
    				for (int j=0;j<vals.length;j++)
    				{
    					Argument a = arg.getArgument(j);
    					if (a._t==ArgumentType.NUMBER_TYPE)
    					{ 
    						NumberArgument narg = (NumberArgument) a; 
    						vals[j] = narg.getNumber().doubleValue(); 
    					}
    				}    				
    				vp.setColorMapValues(vals);
    				vp.repaint();
    			}
    			break;
    			case REDRAW:
    			{
    				int mode = -1;
    				String modeStr = cmd._argv.get(0).toString().toLowerCase(); 
    				if (modeStr.equals("radiate"))
    					mode = RNA.DRAW_MODE_RADIATE;
    				else if (modeStr.equals("circular"))
    					mode = RNA.DRAW_MODE_CIRCULAR;
    				else if (modeStr.equals("naview"))
    					mode = RNA.DRAW_MODE_NAVIEW;
    				else if (modeStr.equals("linear"))
    					mode = RNA.DRAW_MODE_LINEAR;
    				if (mode != -1)
    				  vp.drawRNA(vp.getRNA(), mode);
    			}
    			break;
    			case TOGGLE_SHOW_COLOR_MAP:
    			{
    				vp.setColorMapVisible(!vp.getColorMapVisible());
    			}
    			break;
    			default:
    				throw new Exception(SCRIPT_ERROR_PREFIX+": Method '"+cmd._f+"' unimplemented.");
    		}
    		vp.repaint();
    	}
    }

	
	
	private static Color parseColor(String s)
	{
		Color result = null;
		try {result = Color.decode(s); }
		catch (Exception e) {}
		return result;
	}

	
	private static Vector<Argument> parseArguments(StreamTokenizer st, boolean parType) throws Exception
	{
		Vector<Argument> result = new Vector<Argument>();
		while((st.ttype!=')' && parType) || (st.ttype!=']' && !parType))
		{
			st.nextToken();
			  //System.out.println(""+ (parType?"Par.":"Bra.")+" "+(char)st.ttype);
			switch(st.ttype)
			{
			  case(StreamTokenizer.TT_NUMBER):
			  {
				  result.add(new NumberArgument(st.nval));
			  }
			  break;
			  case(StreamTokenizer.TT_WORD):
			  {
				  Color c = parseColor(st.sval);
				  if (c==null)
					    result.add(new StringArgument(st.sval));
				  else
				    result.add(new ColorArgument(c));
			  }
			  break;
			  case('"'):
			  {
				  result.add(new StringArgument(st.sval));
			  }
			  break;
			  case('['):
			  {
				  result.add(new ArrayArgument(parseArguments(st, false)));
			  }
			  break;
			  case('('):
			  {
				  result.add(new ArrayArgument(parseArguments(st, true)));
			  }
			  break;
			  case(')'):
			  {
				  if (parType)
				    return result;
				  else
					throw new Exception(SCRIPT_ERROR_PREFIX+": Opening "+(parType?"parenthesis":"bracket")+" matched with a closing "+(!parType?"parenthesis":"bracket"));					  
			  }
			  case(']'):
			  {
				  if (!parType)
				    return result;
				  else
					throw new Exception(SCRIPT_ERROR_PREFIX+": Opening "+(parType?"parenthesis":"bracket")+" matched with a closing "+(!parType?"parenthesis":"bracket"));					  
			  }
			  case(','):
				  break;
			  case(StreamTokenizer.TT_EOF):
			  {
				  throw new Exception(SCRIPT_ERROR_PREFIX+": Unmatched opening "+(parType?"parenthesis":"bracket"));
			  }
			  
			}
		}
		return result;
	}
	
	
	private static Command parseCommand(String cmd) throws Exception
	{
		int cut = cmd.indexOf("(");
		if (cut==-1)
		{
			throw new Exception(SCRIPT_ERROR_PREFIX+": Syntax error");
		}
		String fun = cmd.substring(0,cut);
		Function f  = getFunction(fun);
		if (f==Function.UNKNOWN)
		{ throw new Exception(SCRIPT_ERROR_PREFIX+": Unknown function \""+fun+"\""); }
		StreamTokenizer st = new StreamTokenizer(new StringReader(cmd.substring(cut+1)));
		st.eolIsSignificant(false);
		st.parseNumbers();
		st.quoteChar('\"');
		st.ordinaryChar('=');
		st.ordinaryChar(',');
		st.ordinaryChar('[');
		st.ordinaryChar(']');
		st.ordinaryChar('(');
		st.ordinaryChar(')');
		st.wordChars('#', '#');
		Vector<Argument> argv = parseArguments(st,true);
		checkArgs(f,argv);
		Command result = new Command(f,argv); 
		return result;
	}
	
	private static boolean checkArgs(Function f, Vector<Argument> argv) throws Exception
	{
		ArgumentType[] argtypes = getPrototype(f);
		if (argtypes.length!=argv.size())
			throw new Exception(SCRIPT_ERROR_PREFIX+": Wrong number of argument for function \""+f+"\".");
		for (int i=0;i<argtypes.length;i++)
		{
			if (argtypes[i] != argv.get(i)._t)
			{
				throw new Exception(SCRIPT_ERROR_PREFIX+": Bad type ("+argtypes[i]+"!="+argv.get(i)._t+") for argument #"+(i+1)+" in function \""+f+"\".");
			}
		}
		return true;
	}
	
	private static Vector<Command> parseScript(String cmd) throws Exception
	{
		initFunctions();
		Vector<Command> cmds = new Vector<Command>();
		String[] data = cmd.split(";");
		for (int i=0;i<data.length;i++)
		{
			cmds.add(parseCommand(data[i].trim()));
		}
		return cmds;
	}
	
	
}
