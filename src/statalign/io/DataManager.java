package statalign.io;

import java.io.File;
import java.io.IOException;
import java.util.List;

import statalign.base.Utils;
import statalign.io.input.DataReader;

public class DataManager {
	
	List<DataReader> readers;
	
	public void init() {
		readers = Utils.findPlugins(DataReader.class);
	}
	
	/**
	 * Tries to read file using all available input plugins. First tries those that claim to
	 * explicitly support the file extension, that failing tries all other plugins. 
	 * @param file the input file to read
	 * @return the data read or <code>null</code> if all plugins failed to read the file
	 */
	public DataType read(File file) {
		String ext = getExtension(file);
		DataType data;
		for(DataReader reader : readers) {	
			if(reader.supportedExtensions().contains(ext)) {
				try {
					data = reader.read(file);
					System.out.println("Successfully read "+file.getName()+" using "+reader.toString());
					return data;
				} catch (IOException e) {
					System.out.println(e.getMessage());
				}
			}
		}
		for(DataReader reader : readers) {
			if(!reader.supportedExtensions().contains(ext)) {
				try {
					data = reader.read(file);
					System.out.println("Successfully read "+file.getName()+" using "+reader.toString());
					return data;
				} catch (IOException e) {
					System.out.println(e.getMessage());
				}
			}
		}
		return null;
	}
	
	static String getExtension(File file) {
		String name = file.getName();
		int ind = name.lastIndexOf('.');
		if(ind == -1)
			return "";
		return name.substring(ind+1).toLowerCase();
	}
}
