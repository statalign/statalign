package statalign.io.input;

import java.io.IOException;

public class IllegalFormatException extends IOException {
	
	private static final long serialVersionUID = 1L;

	public IllegalFormatException() {
	}
	
	public IllegalFormatException(String str) {
		super(str);
	}
}
