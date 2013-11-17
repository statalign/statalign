package statalign.io.input;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.List;

import statalign.io.DataType;

public abstract class DataReader {

	protected String filename;
	/**
	 * DataReader plugins must override this to return a list of supported file extensions.
	 * @return a list of all lower case file extensions supported by this plugin
	 */
	public abstract List<String> supportedExtensions();
	
	public DataType read(File file) throws IOException {
		filename = file.getName();
		return read(new FileReader(file));
	}

	public DataType read(String fileName) throws IOException {
		filename = fileName;
		return read(new File(fileName));
	}
	
	public abstract DataType read(Reader reader) throws IOException;
	
	@Override
	public String toString() {
		return getClass().getSimpleName();
	}
	
}
