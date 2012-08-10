package statalign.postprocess.utils;

public class EntropyContainer {

	public double obsEntropy;
	public double expEntropy;
	public double sampleEntropy;
	
	
	public EntropyContainer(double obsEntropy, double expEntropy, double sampleEntropy) {
		this.obsEntropy = obsEntropy;
		this.expEntropy = expEntropy;
		this.sampleEntropy = sampleEntropy;
	}
	
	public void setContents(double obsEntropy, double expEntropy) {
		this.obsEntropy = obsEntropy;
		this.expEntropy = expEntropy;
	}
	
}
