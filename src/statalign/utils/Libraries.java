package statalign.utils;

import java.net.URI;
import java.net.URISyntaxException;


public enum Libraries {
	STATALIGN("StatAlign", statalign.StatAlign.version, statalign.StatAlign.webPageURL),
	COLT("Colt", "v1.2.0", "http://acs.lbl.gov/software/colt/"),
	CMATH("Commons Math", "v3.0", "http://commons.apache.org/proper/commons-math/"),
	JAMA("JAMA", "v1.0.3", "http://math.nist.gov/javanumerics/jama/"),
	MPJ("MPJ Express", "v0.38", "http://mpj-express.org/"),
	PPFOLD("PPfold", "v3.0", "http://www.daimi.au.dk/~compbio/pfold/downloads.html"),
	VARNA("VARNA", "v3.8", "http://varna.lri.fr/");
	
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
