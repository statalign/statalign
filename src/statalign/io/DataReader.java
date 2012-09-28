package statalign.io;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.List;

public abstract class DataReader<T> {
	
	public abstract List<String> supportedExtensions();
	
	public T read(File file) throws IOException {
		return read(new FileReader(file));
	}

	public T read(String fileName) throws IOException {
		return read(new File(fileName));
	}
	
	public abstract T read(Reader reader) throws IOException;
	
}
