package statalign.io.input;
//import java.io.BufferedReader;
//import java.io.FileNotFoundException;
//import java.io.FileReader;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.regex.Matcher;
//import java.util.regex.Pattern;
/**
 * Currently unused, this class will be useful in the future console version
 * for reading configuration files.
 */
public class ConfReader {
	//private BufferedReader file;
	//private HashMap<String,String> keyValMap;
	//private ArrayList<String> sectArr;
	/*
	private ConfReader(String fname) throws FileNotFoundException {
		file = new BufferedReader(new FileReader(fname));
		keyValMap = new HashMap<String,String>(30);
		sectArr = new ArrayList<String>(15);
	}
	*/
	/*
	void readConf() throws IOException {
		String line;
		Pattern patSect = Pattern.compile("^\\s*\\[([^\\]]+)\\]\\s*$");
		Pattern patKeyVal = Pattern.compile("^\\s*([\\S&&[^=]]+)\\s*=\\s*\"?(.*?)\"?\\s*$");
		while((line = file.readLine()) != null) {
			String[] sArr = line.split("\\s*[#;]");
			if(sArr.length == 0)
				continue;
			line = sArr[0];
			Matcher m;
			if((m=patSect.matcher(line)).matches())
				sectArr.add(m.group(1));
			else if((m=patKeyVal.matcher(line)).matches())
				keyValMap.put(m.group(1), m.group(2));
		}
	}
	*/
	/*
	void close() throws IOException {
		file.close();
	}
	*/
	/*
	public String getStrVal(String key) {
		return keyValMap.get(key);
	}
	
	public Double getDblVal(String key) {
		String val = keyValMap.get(key);
		return val==null?null:Double.valueOf(val);
	}
	
	public Integer getIntVal(String key) {
		String val = keyValMap.get(key);
		return val==null?null:Integer.valueOf(val);
	}
	
	public String[] getStrArrVal(String key) {
		String val = keyValMap.get(key);
		return val==null?null:val.split("\\s*[,\\s]\\s*");
	}
	
	public double[] getDblArrVal(String key) {
		String[] arr = getStrArrVal(key);
		if(arr == null)
			return null;
		double[] ret = new double[arr.length];
		for(int i = 0; i < arr.length; i++)
			ret[i] = Double.parseDouble(arr[i]);
		return ret;
	}

	public int[] getIntArrVal(String key) {
		String[] arr = getStrArrVal(key);
		if(arr == null)
			return null;
		int[] ret = new int[arr.length];
		for(int i = 0; i < arr.length; i++)
			ret[i] = Integer.parseInt(arr[i]);
		return ret;
	}
	
	public String[] getSectArr() {
		Object[] arr = sectArr.toArray();
		String[] ret = new String[arr.length];
		System.arraycopy(arr, 0, ret, 0, arr.length);
		return ret;
	}
	*/
}
