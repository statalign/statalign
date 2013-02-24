package statalign.utils;

import java.net.URI;
import java.net.URISyntaxException;


public enum Libraries {
	STATALIGN("StatAlign", statalign.StatAlign.version, statalign.StatAlign.webPageURL),
	PPFOLD("PPfold", "v3.0", "http://www.daimi.au.dk/~compbio/pfold/downloads.html"),
	VARNA("VARNA", "v3.8", "http://varna.lri.fr/"),
	MPJ("MPJ Express", "v0.38", "http://mpj-express.org/"),
	JAMA("JAMA", "v1.0.3", "http://math.nist.gov/javanumerics/jama/");
	
	public String longName;
	public String version;
	public String webPage;
	
	private Libraries(String longName, String version, String webPage) {
		this.longName = longName;
		this.version = version;
		this.webPage = webPage;
	}

	public URI getWebPageURI() {
		try {
			return new URI(webPage);
		} catch (URISyntaxException e) {
			return null;
		}
	}
	
	public String getLicenseFile() {
		return "license/"+name()+".txt"; 
	}
	
}
