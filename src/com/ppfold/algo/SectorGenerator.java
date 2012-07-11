package com.ppfold.algo;

/**
 * Generates the sector structure on the submitter machine.
 * 
 * @author Z.Sukosd
 * @see Sector
 */

public class SectorGenerator {

	public static Sector GenerateSectors(int sectorNumber, int distance,
			int divisions, int seqlength, boolean diffbp) {

		Sector sector[];
		sector = new Sector[sectorNumber];

		int rownumber = 0; // refers to job place in triangle
		int columnnumber = 1; // refers to job place in triangle

		int dependencies = 0; // since dependency is always on two adjacent
		// jobs, specify only the first
		int firstSectorInRow = 0; // the ID of the first job in the current row
		// of job
		int lastSectorInRow = divisions; // the ID of the last job in current
		// row of job

		for (int i = 0; i <= sectorNumber - 1; i++) {
			sector[i] = new Sector();
		}

		for (int i = 0; i <= sectorNumber - 1; i++) {
			sector[i].sectorid = i;
			sector[i].dim = distance;
			sector[i].diffbp = diffbp;
			sector[i].initializeBasepairs();
			if(diffbp){
				sector[i].initializeBasepairs2();
			}

			sector[i].fullseqlength = seqlength;
			if (i < sectorNumber - 1) {
				sector[i].next = sector[i + 1];
			} else {
				sector[i].next = null;
			}

			// set dependency
			// we make the triangle "complete"; check later if points are inside
		}
		for (int i = 0; i <= sectorNumber - 1; i++) {
			if (columnnumber > lastSectorInRow) { // if last job in row
				rownumber++;
				columnnumber = 1;
				dependencies = firstSectorInRow;
				firstSectorInRow = i;
				lastSectorInRow = divisions - rownumber;
			}

			if (firstSectorInRow == 0) { // we're in first row
				sector[i].setDependency(null, null);
				int above = i + divisions;
				int jobtoside = i - 1;
				int diagabove = -1;
				if (lastSectorInRow == columnnumber) {
					above = -1;
				}
				if (jobtoside < firstSectorInRow) {
					diagabove = -1;
				} else {
					diagabove = jobtoside + divisions;
				}

				if (above == -1 && diagabove != -1) {
					sector[i].setNextJobs(null, sector[diagabove]);
				} else if (above != -1 && diagabove == -1) {
					sector[i].setNextJobs(sector[above], null);
				} else if (above == -1 && diagabove == -1) {
					sector[i].setNextJobs(null, null);
				} else {
					sector[i].setNextJobs(sector[above], sector[diagabove]);
				}

			} else {
				sector[i].setDependency(sector[dependencies],
						sector[dependencies + 1]);
				dependencies++;
				int above = i + lastSectorInRow;
				int jobtoside = i - 1;
				int diagabove = -1;
				if (columnnumber == lastSectorInRow) {
					above = -1;
				}
				if (jobtoside < firstSectorInRow) {
					diagabove = -1;
				} else {
					diagabove = jobtoside + lastSectorInRow;
				}

				if (above == -1 && diagabove != -1) {
					sector[i].setNextJobs(null, sector[diagabove]);
				} else if (above != -1 && diagabove == -1) {
					sector[i].setNextJobs(sector[above], null);
				} else if (above == -1 && diagabove == -1) {
					sector[i].setNextJobs(null, null);
				} else {
					sector[i].setNextJobs(sector[above], sector[diagabove]);
				}
			}
			// set sector ranges, also preparation for sequence ranges
			int mini = (columnnumber - 1) * distance;
			int minj = rownumber * distance - 1;
			sector[i].setPos(mini, minj);

			if (mini > seqlength - 1) {
				sector[i].seqlength = 0; // set sequence of corresponding job
			} else {
				int seqend = mini + minj + distance;
				if (seqend > seqlength) {
					seqend = seqlength;
				}
				sector[i].seqlength = seqend - mini; // set sequence of
				// corresponding job
			}

			columnnumber++; // take next job in row.

		}

		return sector[sectorNumber - 1];
	}

}
