/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    or obtain one at http://www.gnu.org

    contact:
    address: R.E. Voorrips
            Wageningen University and Research Centre - Plant Breeding
            P.O. Box 16
            6700 AA Wageningen
            The Netherlands
    e-mail: roeland.voorrips@wur.nl
*/

package PedigreeSim;

import JSci.maths.statistics.TDistribution;
import java.io.*;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 *
 * @author Roeland Voorrips
 */
public class Main {

    public static String appName = "PedigreeSim";

    /**
     * @param args the command line arguments
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        simulate(args);
        //Tools.recombtable();
    } //main

    public static void simulate(String[] args) throws Exception {
        File parFile = getParameterFile(args);
        if (parFile==null) { return; }
        System.out.println("Parameter file="+parFile.getName());
        String err = "";
        HashMap<String, String> parMap = readParameterFile(parFile);
        boolean test = false;
        int testIter = 1;
        boolean bivalBidir = false;
        if (parMap.containsKey(strTest)) {
            test = toBoolean(parMap.get(strTest));
            if (parMap.containsKey(strTestIter))
                testIter = Integer.parseInt(parMap.get(strTestIter));
            if (parMap.containsKey(strBivalentsBidirectional))
                bivalBidir = toBoolean(parMap.get(strBivalentsBidirectional));
        }
        boolean interf = false;
        int ploidy;
        if (parMap==null || !parMap.containsKey(strPloidy)) {
            err = parFile.getName()+": no valid "+strPloidy + " found";
        } else
        if (parMap.containsKey(strMapfunction)) {
            interf = parMap.get(strMapfunction).toUpperCase().equals(strKosambi);
        }
        else {
            if(!parMap.containsKey(strHaplostruct)) {
                //mapfunction needed but not present
                err = parFile.getName()+": no valid "+strMapfunction + " found";
            }
        }
        if (!err.isEmpty()) { 
            System.out.println(err);
            return;
        }
        //all data for PopulationData are available
        ploidy = Integer.parseInt(parMap.get(strPloidy));
        double paralMv = 0.0;
        if (parMap.containsKey(strParalMultivalents)) {
            paralMv = Double.parseDouble(parMap.get(strParalMultivalents));
        }
        double pairedCentro = 0.0;
        if (parMap.containsKey(strPairedCentromeres)) {
            pairedCentro = Double.parseDouble(parMap.get(strPairedCentromeres));
        }
        boolean allowNoChiasmata = true;
        if (parMap.containsKey(strAllowNoChiasmata)) {
            allowNoChiasmata = toBoolean(parMap.get(strAllowNoChiasmata));
        }
        boolean naturalPairing = true;
        if (parMap.containsKey(strNaturalPairing)) {
            naturalPairing = toBoolean(parMap.get(strNaturalPairing));
        }
        String missing = strDefaultMissing;
        if (parMap.containsKey(strMissing)) {
            missing = parMap.get(strMissing);
        }
        int seed;
        if (parMap.containsKey(strRandomseed)) {
            seed =  Integer.parseInt(parMap.get(strRandomseed));
        }
        else seed = new Random().nextInt();
        if (!parMap.containsKey(strChromfile)) {
            err = parFile.getName()+": no valid "+strChromfile + " found";
        }
        if (!parMap.containsKey(strOutput)) {
            err = parFile.getName()+": no valid "+strOutput + " found";
        } else {
            if (!canCreateFile(parMap.get(strOutput)+".t64")) {
                err = parFile.getName()+": cannot create "+strOutput + ".*";
            }
        }
        if (!err.isEmpty()) {
            System.out.println(err);
            return;
        }
        PopulationData popdata = new PopulationData(ploidy, interf,
                paralMv, pairedCentro, naturalPairing, allowNoChiasmata,
                bivalBidir, seed, missing);
        
        try { 
            readChromosomeFile(parMap.get(strChromfile), popdata); 
            if (popdata.chromCount()==0)
                throw new Exception(strChromfile+": no chromosomes");
            if (!popdata.allowNoChiasmata) {
                double minLength = 0.70; //minimum chromosome length if allowNoChiasmata is true
                for (Chromosome chr: popdata.getChromosome()) {
                    if (chr.getLength()<minLength)
                        throw new Exception("if " + strAllowNoChiasmata +
                                " is false the minimum chromosome length is " +
                                (minLength*100.0) + " cM");
                }
            }
            if (parMap.containsKey(str2nGametes)) {
                //cannot be set in the constructor (for safety), 
                //therefore added after construction of popdata
              if (parMap.get(str2nGametes).toUpperCase().equals(strFDR)) {
                  popdata.unreducedGametes = 1; //FDR
              } else popdata.unreducedGametes = 2; //SDR      
            }
        }    
        catch (Exception ex) {
            System.out.println(ex.getMessage());
            return;
        }

        //read or generate&write pedigree:

        ArrayList<String[]> pedList;

        if (test) {
            popdata.testMode = true;
            if (parMap.containsKey(strPrintGametes)) {
                popdata.testPrintGametes = Integer.parseInt(parMap.get(strPrintGametes));
            }
            if (parMap.containsKey(strPrintEachIter)) {
                popdata.testPrintEachIter = toBoolean(parMap.get(strPrintEachIter));
            }
            if (parMap.containsKey(strPrintMapdistances)) {
                popdata.testPrintMapdistances = toBoolean(parMap.get(strPrintMapdistances));
            }
            if (parMap.containsKey(strPrintPubTables)) {
                popdata.testPrintPubTables = toBoolean(parMap.get(strPrintPubTables));
            }
            //create population of just one individual:
            String[] TestParent = new String[] {"Parent", missing, missing};
            pedList = new ArrayList<String[]>();
            pedList.add(TestParent);
            err = popdata.createPopulation(pedList);
            if (err.isEmpty()) popdata.getIndiv(0).setHaploStruct(
                    Individual.calcFounderHaplostruct(-1, popdata));
            if (parMap.containsKey(strPopsize)) {
                popdata.testMeioseCount =
                        Integer.parseInt(parMap.get(strPopsize));
            }
            else err = strPopsize + " not found in parameter file";
        }
        else {
            if (parMap.containsKey(strPedfile)) {
                //get individuals from pedigree file
                ArrayList<String[]> pedListOrig = readPedigreeFile(parMap.get(strPedfile));
                if (!pedListOrig.get(0)[0].isEmpty()) {
                    System.out.println(""+pedListOrig.get(0)[0]);
                    return;
                }
                pedListOrig.remove(0);
                try {
                    pedList = popdata.sortPedigree(pedListOrig);
                    err = popdata.createPopulation(pedList);
                    //check if any changes were made:
                    int i = pedList.size()-1;
                    while (i>=0 && pedList.get(i)[0].equals(pedListOrig.get(i)[0])) {
                        i--;
                    }
                    if (i>=0) {
                        //changes due to sorting, write the new pedigree to file:
                        writePedigreeFile(parMap.get(strOutput)+".ped",popdata);
                    }
                } catch (Exception ex) {
                    System.out.println("Error in pedigree: "+ex.getMessage());
                    return;
                }
            }
            else {
                if (!parMap.containsKey(strPoptype) ||
                    !parMap.containsKey(strPopsize) ) {
                    err = "Either a valid "+strPedfile+" or both "+strPoptype+" and "+
                            strPopsize+" must be specified";
                }
                else {
                    pedList = generateStandardPopulation(
                            parMap.get(strPoptype).toUpperCase(),
                            Integer.parseInt(parMap.get(strPopsize)),
                            missing);
                    if (pedList.isEmpty()) {
                        err = "Invalid population type or population size==0";
                    }
                    else {
                        err = popdata.createPopulation(pedList);
                        if (err.isEmpty())
                            writePedigreeFile(parMap.get(strOutput)+".ped",popdata);
                    }
                }
            } //if PEDFILE else
        } // if test else

        /* pedigree is now available;
         * read or simulate the HaploStruct
         * and generate output based on map and founder genotypes
         */

        if (err.isEmpty()) {

            if (test) {
                if (parMap.containsKey(strHaplostruct)) {
                    try {
                        readHaploStructFiles(parMap.get(strHaplostruct),popdata);
                    } catch (Exception ex) {
                        System.out.println(ex.getMessage());
                        return;
                    }
                }
                // else
                if (!parMap.containsKey(strMapfile)) {
                    err = strMapfile + " not found in parameter file";
                }
                else {
                    try {
                        readMapFile(parMap.get(strMapfile),popdata);
                    } catch (Exception ex) {
                        err = "Error reading "+parMap.get(strMapfile)+": "+ex.getMessage();
                    }
                }
                if (err.isEmpty()) {
                    popdata.testIter = testIter;
                    testMeiosis(popdata,parMap.get(strOutput),2.63);
                }
            } //if test
            else { //normal simulation

                //read and/or simulate&write haplostructs:
                int maxfounderallele = -1;

                if (parMap.containsKey(strHaplostruct)) {
                    try {
                        maxfounderallele =
                                readHaploStructFiles(parMap.get(strHaplostruct),
                                                     popdata);
                    } catch (Exception ex) {
                        System.out.println(ex.getMessage());
                        return;
                    }
                }
                /* three cases:
                 * some individuals in haplostructs but not in pedigree: error
                 * some individuals in pedigree but not in haplostructs: simulate
                 *      new haplostructs for those (also in the standard case
                 *      where no haplostructs files provided: then all pedigree
                 *      individuals are simulated)
                 * individuals in haplostructs match those in pedigree exactly:
                 *      no simulation of haplostructs needed
                 * In both last cases actual genotypes are generated if
                 * map and founder genotypes available
                 */
                //simulate new haplostruct for pedigree
                simulateHaploStructs(popdata, maxfounderallele);
                writeHaploStructFiles(parMap.get(strOutput),popdata);
                if (!err.isEmpty()) {
                    System.out.println(err);
                    return;
                }

                //do we have to produce genotypes?
                if (parMap.containsKey(strMapfile)) {
                    //produce and output genotypes
                    try {
                        readMapFile(parMap.get(strMapfile),popdata);
                    } catch (Exception ex) {
                        err = "Error reading "+parMap.get(strMapfile)+": "+ex.getMessage();
                    }
                    boolean founderGenotypes = false;
                    if (err.isEmpty()) {
                        if (parMap.containsKey(strFounderfile)) {
                            try {
                                readFounderGenotypesFile(parMap.get(strFounderfile),popdata);
                                founderGenotypes = true;
                            } catch (Exception ex) {
                                err = "Error reading "+parMap.get(strFounderfile)+": "+ex.getMessage();
                            }
                        }
                    }
                    if (err.isEmpty()) {
                        try {
                            int whichAllele;
                            if (founderGenotypes) {
                                //actual alleles present, write file with actual allele genotypes:
                                writeGenotypesFile(parMap.get(strOutput)+"_genotypes.dat",false,popdata);
                                whichAllele = MAX_ALLELE;
                            }
                            else whichAllele = 0;
                            //write file with founderallele genotypes:
                            writeGenotypesFile(parMap.get(strOutput)+"_founderalleles.dat",true,popdata);
                            //write allele dosage file for either the MAX_ALLELE is actual alleles are present,
                            //and else for founderallele 0:
                            writeAlleledoseFile(parMap.get(strOutput)+"_alleledose.dat",whichAllele,founderGenotypes,popdata);
                            //write allAlleledose.dat:
                            //- for all actual alleles only if some or all loci have more than 2 alleles
                            //- or for all founder alleles
                            if (!founderGenotypes || popdata.getMaxAlleleCount()>2) {
                                writeAllAlleledoseFile(parMap.get(strOutput)+"_allAlleledose.dat",founderGenotypes,popdata);
                            }    
                        } catch (Exception ex) {
                            err = "Error writing files: "+ex.getMessage();
                        }
                    }
                }
                else {
                    //output possible error:
                    //Haplostruct files given but no Map: nothing to do
                    if ( parMap.containsKey(strHaplostruct) &&
                         !parMap.containsKey(strMapfile) ) {
                        err = parFile.getName()+": when "+strHaplostruct+" is given also "+
                                strMapfile+" must be specified";
                    }
                }
            } // if test else
        }
        
        if (!err.isEmpty()) {
            System.out.println(err);
        }

    } //simulate

    private static final char commentChar = ';';
    //Standard keys:
    private static final String strMissing = "MISSING";
    private static final String strRandomseed = "SEED";
    private static final String strPloidy = "PLOIDY";
    private static final String strMapfunction = "MAPFUNCTION";
    private static final String strChromfile = "CHROMFILE";
    private static final String strPedfile = "PEDFILE";
    private static final String strMapfile = "MAPFILE";
    private static final String strFounderfile = "FOUNDERFILE";
    private static final String strHaplostruct = "HAPLOSTRUCT";
    private static final String strOutput = "OUTPUT";
    private static final String strPoptype = "POPTYPE";
    private static final String strPopsize = "POPSIZE";
    //Advanced keys:
    private static final String strParalMultivalents = "PARALLELMULTIVALENTS";
    private static final String strPairedCentromeres = "PAIREDCENTROMERES";
    private static final String strAllowNoChiasmata = "ALLOWNOCHIASMATA";
    private static final String strNaturalPairing = "NATURALPAIRING";
    //Test mode keys:
    private static final String strTest = "TEST";
    private static final String strTestIter = "TESTITER";
    private static final String strBivalentsBidirectional = "BIVALENTSBIDIRECTIONAL";
    private static final String strPrintGametes = "PRINTGAMETES";
    private static final String strPrintEachIter = "PRINTEACHITER";
    private static final String strPrintMapdistances = "PRINTMAPDISTANCES";
    private static final String strPrintPubTables = "PRINTPUBTABLES"; //only in combination with specific chrom and map files 
    private static final String str2nGametes = "2NGAMETES";
    // Value strings:
    private static final String strDefaultMissing = "NA";
    private static final String strHaldane = "HALDANE";
    private static final String strKosambi = "KOSAMBI";
    private static final String strF1 = "F1"; //F1 population from 2 heterozygous parents
    private static final String strF2 = "F2"; //F2 population from 2 heterozygous perents
    private static final String strBC = "BC"; //BC population from 2 heterozygous parents
    private static final String strS1 = "S1"; /*progeny from selfing one heterozygous parent
                                         *(equivalent to F2 from two homozygous parents
                                         * for observed alleles, but not for founder
                                         * alleles)
                                         */
    private static final String strFDR = "FDR";
    private static final String strSDR = "SDR";
            
    
    public static HashMap<String, String> readParameterFile(File parFile)  {
        if (parFile==null) return null;
        HashMap<String, String> parMap = new HashMap<String, String>();
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(parFile));
        } catch (Exception ex) { return null; }
        do {
            String s;
            try {
                s = in.readLine();
            } catch (Exception ex) { s = null; }
            if (s==null) break;
            String[] keyval = getKeyandValue(s);
            if ( ( keyval[0].equals(strPloidy) &&
                    isValidPloidy(keyval[1]) )
                    ||
                 ( keyval[0].equals(strMapfunction) &&
                    (keyval[1].toUpperCase().equals(strHaldane) ||
                     keyval[1].toUpperCase().equals(strKosambi)) )
                     ||
                 ( (keyval[0].equals(strChromfile) ||
                    keyval[0].equals(strPedfile) ||
                    keyval[0].equals(strMapfile) ||
                    keyval[0].equals(strFounderfile)) &&
                      new File(keyval[1]).canRead() )
                     ||
                 ( keyval[0].equals(strHaplostruct) &&
                    (new File(keyval[1]+".hsa").canRead() &&
                     new File(keyval[1]+".hsb").canRead()) )
                     ||
                 ( keyval[0].equals(strOutput) &&
                    !keyval[1].isEmpty() ) //further check needed when (re)writing the files
                     ||
                 ( keyval[0].equals(strPoptype) &&
                    (keyval[1].equals(strF1) ||
                     keyval[1].equals(strF2) ||
                     keyval[1].equals(strBC) ||
                     keyval[1].equals(strS1)) ) //correction of Fabian Grandke
                     ||
                 ( (keyval[0].equals(strPopsize) ||
                    keyval[0].equals(strTestIter))&&
                    isPositiveInteger(keyval[1]))
                     ||
                 ( (keyval[0].equals(strAllowNoChiasmata) ||
                    keyval[0].equals(strNaturalPairing) ||
                    keyval[0].equals(strBivalentsBidirectional) ||
                    keyval[0].equals(strTest) ||
                    keyval[0].equals(strPrintEachIter) ||
                    keyval[0].equals(strPrintMapdistances) ||
                    keyval[0].equals(strPrintPubTables) ) &&
                     isBoolean(keyval[1]) )
                     ||
                 ( (keyval[0].equals(strParalMultivalents) ||
                    keyval[0].equals(strPairedCentromeres)) &&
                    isProbability(keyval[1]) )
                     ||   
                 ( keyval[0].equals(strMissing) &&
                    !keyval[1].isEmpty() )
                     ||
                 ( (keyval[0].equals(strRandomseed) ||
                    keyval[0].equals(strPrintGametes) ) &&
                    isInteger(keyval[1]) ) 
                    ||
                 ( keyval[0].equals(str2nGametes) &&
                    (keyval[1].toUpperCase().equals(strFDR) ||
                     keyval[1].toUpperCase().equals(strSDR)) )
               ) {
                        parMap.put(keyval[0], keyval[1]);
            }
        } while(true);
        return parMap;
    } //readParameterFile

    public static void readChromosomeFile(String fName, PopulationData popdata)
            throws FileNotFoundException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 || !words[0].toUpperCase().equals("CHROMOSOME") ||
                !words[1].toUpperCase().equals("LENGTH") ||
                !words[2].toUpperCase().equals("CENTROMERE")) {
            throw new Exception(fName+": header line not found");
        }
        //Read data lines:
        int wordsNeeded = popdata.ploidy==2 ? 3 : 5;
        do {
            s=in.readLine(); //skip headerline
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<wordsNeeded) {
                    throw new Exception(fName+": insufficient data on line '"+s+"'");
                } 
                int i = wordsNeeded-1;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception(fName+": missing data are not allowed");
                }
                String chromName = words[0];
                double length = Double.parseDouble(words[1])/100.0; //cM to Morgan;
                double centromerePos = Double.parseDouble(words[2])/100.0; //cM to Morgan;
                if (popdata.ploidy==2) {
                    popdata.addChromosome(new Chromosome(chromName,length,
                            centromerePos,popdata));
                }
                else {
                    double prefPairingProb = Double.parseDouble(words[3]);
                    double fracQuadrivalents = Double.parseDouble(words[4]);
                    popdata.addChromosome(new TetraploidChromosome(chromName,length,
                            centromerePos,prefPairingProb,fracQuadrivalents,
                            popdata));
                }
            }
        } while(true);
    } //readChromosomeFile

    public static void readMapFile(String fName, PopulationData popdata)
            throws FileNotFoundException, IOException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 || !words[0].toUpperCase().equals("MARKER") ||
                !words[1].toUpperCase().equals("CHROMOSOME") ||
                !words[2].toUpperCase().equals("POSITION")) {
            throw new Exception(fName+": header line not found");
        }
        //Read data lines:
        String lastName="";
        Chromosome chrom = null;
        String[] alleleNames = new String[] {"0","1"};
        Locus locus;
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<3) {
                    throw new Exception("Too few items on line '"+s+"'");
                }
                int i = 2;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception(fName+": missing data are not allowed");
                }
                String mrkName = words[0];
                String chromName = words[1];
                double pos = Double.parseDouble(words[2])/100.0; //cM to Morgan
                if (!chromName.equals(lastName)) {
                    int chr = 0;
                    do {
                        chrom = popdata.getChrom(chr);
                        chr++;
                    } while(chr<popdata.chromCount() &&
                            !chromName.equals(chrom.getChromName()));
                    if (!chromName.equals(chrom.getChromName())) {
                        throw new Exception("chromName "+chromName+" not present");
                    }
                    lastName = chromName;
                }
                locus = new Locus(mrkName,pos,alleleNames,popdata);
                chrom.addLocus(locus);
            }
        } while(true);
    } //readMapFile

    public static void writePedigreeFile(String fName, PopulationData popdata) throws IOException {
        PrintWriter out = new PrintWriter(
                new BufferedWriter(new FileWriter(fName)));
        out.println("Name\tParent1\tParent2");
        for (Individual ind: popdata.getIndividual()) {
            out.print(ind.getIndivName());
            if (ind.isFounder()) {
                out.println("\t"+popdata.missing+"\t"+popdata.missing);
            }
            else {
                out.println("\t"+ind.getParents()[0].getIndivName()+
                            "\t"+ind.getParents()[1].getIndivName());
            }
        }
        out.flush();
    } //writePedigreeFile

    /**
     * readPedigreeFile
     * @param fName
     * @return an ArrayList<String[]> of which the first item is an error
     * @throws IOException Only if readline doesn't work, else error message
     * message ("" if no error) and all following items have 3 Strings:
     * individual name and names of parent1 and parent 2 (both parents should be
     * strMissing if individual is founder)
     */
    public static ArrayList<String[]> readPedigreeFile(String fName)
             throws IOException {
        ArrayList<String[]> result = new ArrayList<String[]>();
        result.add(new String[] {""}); //first item is error message
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            result.get(0)[0] =fName+": cannot open file";
            return result;
        }
        String s;
        //find the header line:
        do {
            s = in.readLine();
            if (s==null) break;
        } while (s.trim().isEmpty() ||
                 s.trim().charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 ||
            !words[0].toUpperCase().equals("NAME") ||
            !words[1].toUpperCase().equals("PARENT1") ||
            !words[2].toUpperCase().equals("PARENT2") ) {
            result.get(0)[0] = fName+": header line not found";
            return result;
        }
        //read the first 3 words of each pedigree line
        do {
            s = in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<3 ||
                    words[1].charAt(0)==commentChar ||
                    words[2].charAt(0)==commentChar) {
                    result.get(0)[0] = "Incomplete line in "+fName+": "+s;
                    return result;
                }
                result.add(new String[] {words[0], words[1],words[2]});
            }
        } while (true);
        return result;
    } //readPedigreeFile

    public static final char HOMOLOG_SEPARATOR = '_';

    private static boolean testGenotypesfileCaption(String s, 
            int hspos, Individual ind, int homolog) {
        if (s==null || s.length()<3) return false;
        if (s.charAt(hspos) != HOMOLOG_SEPARATOR ||
                !s.substring(0, hspos).equals(ind.getIndivName()))
            return false;
        try {
            int p = Integer.parseInt(s.substring(hspos + 1));
            return p == homolog;
        } catch (Exception ex) {return false;}
    } //testGenotypesfileCaption

    public static void readFounderGenotypesFile(String fName, PopulationData popdata)
            throws FileNotFoundException, IOException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length==0 || !words[0].toUpperCase().equals("MARKER")) {
            throw new Exception("header line not found");
        }
        //read headerline:
        //founders may appear in different order than in popdata.individuals
        if (words.length != 1 + popdata.founderAlleleCount) {
            throw new Exception("must contain one column for the marker name and ploidy columns for each founder");
        }
        if (!words[0].toUpperCase().equals("MARKER")) {
            throw new Exception("Header doesn't start with 'marker'");
        }
        int nf = popdata.founderAlleleCount / popdata.ploidy; //number of founders
        int[] indnr = new int[nf]; //will contain the individual number for each
                                   //founder in order of header
        for (int i=0; i<nf; i++) {
            String iname = words[1+i*popdata.ploidy];
            int hspos = iname.lastIndexOf(HOMOLOG_SEPARATOR);
            Individual ind = 
                    popdata.getIndiv(iname.substring(0, hspos));
            if (ind==null || ind.getParents()[0]!=null) {
                throw new Exception(iname.substring(0, hspos)+" is not a founder");
            }
            indnr[i] = popdata.getIndividual().indexOf(ind);
            for (int hom=0; hom<popdata.ploidy; hom++) { 
                iname = words[1+i*popdata.ploidy+hom];
                if (!testGenotypesfileCaption(iname, hspos, ind, hom+1)) {
                    throw new Exception(iname+" found but "+
                            ind.getIndivName()+HOMOLOG_SEPARATOR+(hom+1)+
                            " expected");
                }
            }
        }
        // now all captions seem ok
        // if no duplicates occur all founders are represented:
        for (int i=1; i<nf; i++) {
            for (int j=0; j<i; j++) {
                if (indnr[j]==indnr[i]) {
                    throw new Exception("Founder "+popdata.getIndiv(indnr[i])+
                            " occurs more than once in header");
                }
            }
        }

        //read the lines per marker, in order of map:
        int cix = 0;
        Chromosome chrom = popdata.getChrom(cix);
        int loc=0;
        int wordsNeeded = 1 + popdata.founderAlleleCount;
        do {
            s = in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length != wordsNeeded) {
                   throw new Exception("not all lines have "+
                           wordsNeeded+" items");
                }
                int i = wordsNeeded -1;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception("missing data are not allowed");
                }
                String mrkName = words[0];
                if (chrom==null ||
                    loc>=chrom.getLocus().size() ||
                    !chrom.getLocus().get(loc).getLocusName().equals(mrkName)) {
                    throw new Exception("invalid marker: "+words[0]);
                }
                String[] alleleNames = new String[popdata.founderAlleleCount];
                for (i=0; i<nf; i++) {
                    int headAllele = popdata.getIndiv(indnr[i]).
                            getHaploStruct(cix,0).
                            getFounderAt(chrom.getHeadPos());
                    for (int hom=0; hom<popdata.ploidy; hom++) {
                        alleleNames[headAllele+hom] = words[1+i*popdata.ploidy+hom];
                    }
                }
                chrom.getLocus().get(loc).setAlleleName(alleleNames);
                loc++;
                if (loc==chrom.getLocus().size()) {
                    loc = 0;
                    cix++;
                    if (cix==popdata.chromCount()) chrom=null;
                    else chrom = popdata.getChrom(cix);
                }
            }
        } while(true);
        if (chrom!=null) {
            throw new Exception(fName+": end of file found while reading chromosome "+
                    chrom.getChromName());
        }
    } //readFounderGenotypesFile

    /**
     * writeGenotypesFile writes a tab-separated text file with for each 
     * individual and each locus on the map the popdata.ploidy alleles.
     * @param fName the file name; an existing file will be overwritten
     * @param founders if true the founder alleles are written, else the
     * corresponding allele names
     * @param popdata
     * @throws IOException
     * @throws Exception 
     */
    public static void writeGenotypesFile(String fName, boolean founders, PopulationData popdata)
            throws IOException, Exception {
        System.out.println("write file "+fName);
        //File outFile = new File(fName);
        //PrintWriter out = new PrintWriter(new FileWriter(outFile));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName)));
        // write header line:
        out.print("marker");
        for (int i=0; i<popdata.indivCount(); i++) {
            for (int h=0; h<popdata.ploidy; h++) {
                out.print("\t"+ popdata.getIndiv(i).getIndivName()+"_"+(h+1));
            }
        }
        out.println();
        for (int c=0; c<popdata.chromCount(); c++) {
            Chromosome chrom = popdata.getChrom(c);
            System.out.println("chrom "+chrom.getChromName());
            ArrayList<Locus> loci = chrom.getLocus();
            for (int loc=0; loc<loci.size(); loc++) {
                //instead of building a line as a string in memeory and printing
                //it at once it is much faster to print every element separately:
                out.print(loci.get(loc).getLocusName());
                for (int i=0; i<popdata.indivCount(); i++) {
                    if (founders) {
                        for (int h=0; h<popdata.ploidy; h++) {
                            out.print("\t"+ popdata.getIndiv(i).getHaploStruct(c,h).
                                    getFounderAt(loci.get(loc).getPosition()));
                        }
                    } else { //not founders but locus alleles
                        String[] indall = popdata.getIndiv(i).getLocusAllele(c, loc);
                        for (int h=0; h<popdata.ploidy; h++) {
                            out.print("\t"+ indall[h]);
                        }
                    }
                }
                out.println();
            }
        }
        out.flush();
    } //writeGenotypesFile

    public static final int MIN_ALLELE=-1;
    public static final int MAX_ALLELE=-2;

    /**
     * writeAlleledoseFile:
     * writes a tab-separated text file with for each individual and each
     * locus on the map the dosage (number of copies) of the allele specified
     * in whichAllele
     * @param fName the file name; an existing file will be overwritten
     * @param whichAllele specifies for which of the alleles the dosage will
     * be written. If actualAlleles is true and the value of whichAllele
     * is MIN_ALLELE or MAX_ALLELE the allele with
     * the lowest or highest allele name will be counted, else the allele
     * corresponding with founder allele whichAllele will be counted (if
     * whichAllele not in 0..popdata.founderAlleleCount-1 the allele
     * corresponding with founderAllele 0 will be counted). 
     * If actualAlleles
     * is false the founder allele corresponding to whichAllele is counted
     * (and if whichAllele is outside 0 ... popdata.founderAlleleCount-1
     * the counts will all be 0)
     * @param actualAlleles if true the counting is based on the actual
     * alleles corresponding to each founder allele; if false it is based
     * on the founder alleles directly (and then whichAllele specifies
     * the founder allele counted)
     * @param popdata
     * @throws IOException
     * @throws Exception 
     */
    public static void writeAlleledoseFile(String fName, int whichAllele,
            boolean actualAlleles, PopulationData popdata)
            throws IOException, Exception {
        System.out.println("write file "+fName);
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName)));
        // write header line:
        out.print("marker");
        for (int i=0; i<popdata.indivCount(); i++) {
            out.print("\t"+ popdata.getIndiv(i).getIndivName());
        }
        out.println();
        for (int c=0; c<popdata.chromCount(); c++) {
            Chromosome chrom = popdata.getChrom(c);
            System.out.println("chrom "+chrom.getChromName());
            ArrayList<Locus> loci = chrom.getLocus();
            for (int loc=0; loc<loci.size(); loc++) {
                out.print(loci.get(loc).getLocusName());
                double locpos = loci.get(loc).getPosition();
                String matchName = "";
                if (actualAlleles) {
                    switch (whichAllele) {
                        case MIN_ALLELE:
                             matchName = loci.get(loc).getMinAlleleName(); break;
                        case MAX_ALLELE:
                             matchName = loci.get(loc).getMaxAlleleName(); break;
                        default: {
                            if (whichAllele<0 || 
                                    whichAllele>=popdata.founderAlleleCount) {
                                whichAllele = 0;
                            }
                            matchName = loci.get(loc).getAlleleName(whichAllele); 
                        } break;
                    }    
                }
                for (int i=0; i<popdata.indivCount(); i++) {
                    int dose = 0;
                    if (actualAlleles) {
                        String[] indall = popdata.getIndiv(i).getLocusAllele(c, loc);
                        for (int h=0; h<popdata.ploidy; h++) {
                            if (matchName.equals(indall[h])) dose++;
                        }
                    }
                    else { // count founder allele whichAllele instead of actual allele
                        int[] indall = popdata.getIndiv(i).getFounderAlleles(c, locpos);
                        for (int h=0; h<popdata.ploidy; h++) {
                            if (whichAllele==indall[h]) dose++;
                        }
                    }
                    out.print("\t" + dose);
                }
                out.println();
            }
        }
        out.flush();
    } //writeAlleledoseFile

    /**
     * writeAllAlleledoseFile:
     * writes a tab-separated text file with for each individual and each
     * locus on the map the dosage (number of copies) of the allele specified
     * in whichAllele
     * @param fName the file name; an existing file will be overwritten
     * @param actualAlleles if true the counting is based on the actual
     * alleles corresponding to each founder allele; if false it is based
     * on the founder alleles directly 
     * @param popdata
     * @throws IOException
     * @throws Exception 
     */
    public static void writeAllAlleledoseFile(String fName, 
            boolean actualAlleles, PopulationData popdata)
            throws IOException, Exception {
        System.out.println("write file "+fName);
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName)));
        // write header line:
        out.print("marker\tallele");
        for (int i=0; i<popdata.indivCount(); i++) {
            out.print("\t"+ popdata.getIndiv(i).getIndivName());
        }
        out.println();
        for (int c=0; c<popdata.chromCount(); c++) {
            Chromosome chrom = popdata.getChrom(c);
            System.out.println("chrom "+chrom.getChromName());
            ArrayList<Locus> loci = chrom.getLocus();
            for (int loc=0; loc<loci.size(); loc++) {
                double locpos = loci.get(loc).getPosition();
                String matchName;
                if (actualAlleles) {
                    String[] alleles = loci.get(loc).getUniqueAlleleNames();
                    for (String allele : alleles) {
                        matchName = allele;
                        out.print(loci.get(loc).getLocusName()+"\t"+matchName);  
                        for (int i=0; i<popdata.indivCount(); i++) {
                            int dose = 0;
                            String[] indall = popdata.getIndiv(i).getLocusAllele(c, loc);
                            for (int h=0; h<popdata.ploidy; h++) {
                                if (matchName.equals(indall[h])) dose++;
                            }
                            out.print("\t" + dose);
                        }
                        out.println();
                    }
                }
                else {
                    //not actualAlleles, count founder alleles
                    for (int a=0; a<popdata.founderAlleleCount; a++) {
                        out.print(loci.get(loc).getLocusName()+"\t"+a);  
                        for (int i=0; i<popdata.indivCount(); i++) {
                            int dose = 0;
                            int[] indall = popdata.getIndiv(i).getFounderAlleles(c, locpos);
                            for (int h=0; h<popdata.ploidy; h++) {
                                if (a==indall[h]) dose++;
                            }
                            out.print("\t" + dose);
                        }
                        out.println();
                    }
                }
            } //for loc
        } // for c
        out.flush();
    } //writeAllAlleledoseFile

    /**
     * 2 files store the HaploStructs of the whole population:
     * *.hsa (the founder segments: set of integers) and
     * *.hsb (the recombination points: set of doubles)
     * for each chromosome the max. length of each haplostruct is used
     * for the whole file
     * @param fName
     * @param popdata
     * @throws IOException
     */
    public static void writeHaploStructFiles(String fName, PopulationData popdata)
            throws IOException {
        PrintWriter hsa = new PrintWriter(
                new BufferedWriter(new FileWriter(fName+".hsa")));
        PrintWriter hsb = new PrintWriter(
                new BufferedWriter(new FileWriter(fName+".hsb")));
        int segcount = 0;
        for (Chromosome chrom: popdata.getChromosome()) {
            int n = chrom.maxSegCount();
            if (n>segcount) segcount = n;
        }

        for (Individual ind: popdata.getIndividual()) {
            for (Chromosome chrom: popdata.getChromosome()) {
                int chrNum = chrom.getChromNumber();
                for (int p=0; p<popdata.ploidy; p++) {
                    hsa.print(ind.getIndivName()+"\t"+chrom.getChromName()+"\t"+(p+1));
                    HaploStruct hs = ind.getHaploStruct(chrNum,p);
                    int hslen = hs.segmentCount();
                    for (int s=0; s<hslen; s++) {
                        hsa.print("\t"+hs.getFounder().get(s));
                        if (s>0) {
                            hsb.print("\t" + hs.getRecombPos().get(s)*100.0); //Morgan to cM
                        }
                    }
                    for (int s=hslen; s<segcount; s++) {
                        hsa.print("\t"+popdata.missing);
                        hsb.print("\t"+popdata.missing);
                    }
                    hsa.println();
                    hsb.println();
                }
            }
        }
        hsa.flush();
        hsb.flush();
    } //writeHaploStructFiles

    public static int readHaploStructFiles(String fName, PopulationData popdata)
            throws Exception {
        //advance information: chromosomes(name, length etc; ploidy; pedigree
        // hsa has the info on individual and chromosome, and the sequences of founder alleles:
        BufferedReader hsa = new BufferedReader(new FileReader(new File(fName+".hsa")));
        // hsb has the sequences of recombination positions:
        BufferedReader hsb = new BufferedReader(new FileReader(new File(fName+".hsb")));
        //first we make sure that all current haplotypes of individuals are removed
        //to allow a check for individuals occurring multiple times in the haplofiles
        //(Note that founders already get haplotypes in their constructor)
        for (Individual ind: popdata.getIndividual()) {
            ind.clearAllHaplostruct();
        }
        String s;
        int line = 0;
        int maxfounderallele = -1;
        String iName = "";
        HaploStruct[][] hs = null;
        do {
            s=hsa.readLine();
            if (s!=null) s = s.trim();
            if (s==null || s.equals("")) break;
            String[] words = Tools.readWords(s);
            if (words.length<4) {
                throw new Exception(fName+".hsa: line "+(line+1)+" has too few elements");
            }
            if (line %(popdata.chromCount()*popdata.ploidy) == 0) {
                Individual ind = popdata.getIndiv(words[0]);
                if (ind==null) {
                    throw new Exception(fName+".hsa: individual "+words[0]+" doesn't exist");
                }
                if (ind.getAllHaploStruct()!=null) {
                    throw new Exception(fName+".hsa: individual "+words[0]+" occurs more than once");
                }
                iName = words[0];
                hs = new HaploStruct[popdata.chromCount()][popdata.ploidy];
            }
            else if (!words[0].equals(iName)) {
                throw new Exception(fName+".hsa/hsb: Insufficient lines for individual "+iName);
            }
            int indline = line % (popdata.chromCount()*popdata.ploidy); //nth line for this individual
            int chrnum = indline / popdata.ploidy;                      //chromosome number
            int homol = indline % popdata.ploidy;                       //homolog number
            Chromosome chrom = popdata.getChrom(chrnum);
            if (!words[1].equals(popdata.getChrom(chrnum).getChromName())) {
                throw new Exception(fName+".hsa/hsb, line "+(line+1)+": chromosome "+
                        chrom.getChromName()+" expected but "+
                        words[1]+" found");
            }
            else {
                int homnr;
                try { homnr = Integer.valueOf(words[2]);
                } catch (Exception ex) {
                    throw new Exception(fName+".hsa/hsb, line "+(line+1)+": homolog "+
                        (homol+1) + " expected, but "+words[2]+" found");
                }
                if (homnr != (homol+1)) {
                    throw new Exception(fName+".hsa/hsb, line "+(line+1)+": homolog "+
                            (homol+1) + " expected, but "+words[2]+" found");
                }
            }
            String hsbline = hsb.readLine();
            if (hsbline==null) {
                throw new Exception(fName+".hsb: less lines than "+fName+".hsa");
            }
            int founderAllele = -1;
            try {
                founderAllele = Integer.valueOf(words[3]);
            } catch (Exception ex) {}
            if (founderAllele<0 || founderAllele>=popdata.founderAlleleCount ) {
                throw new Exception(fName+".hsa, line "+(line+1)+
                        ": invalid founder allele '"+words[3]+"'");
            }
            if (founderAllele > maxfounderallele) maxfounderallele = founderAllele;
            String[] hsbwords = Tools.readWords(hsbline);
            if (hsbwords.length != words.length-4) {
                throw new Exception(fName+".hsa/hsb, line "+(line+1)+
                        ": number of items in both files don't match");
            }
            HaploStruct hast = new HaploStruct(chrom, founderAllele);
            //repeat reading pairs of values from hsa and hsb:
            boolean done = false;
            int r = 0;
            while (!done && r<hsbwords.length) {
                if (words[4+r].equals(popdata.missing)) {
                    done = true;
                    if (!hsbwords[r].equals(popdata.missing)) {
                        throw new Exception(fName+".hsa/hsb, line "+(line+1)+
                                ": number of non-missing items in both files don't match");
                    }
                }
                else {
                    founderAllele = -1;
                    try {
                        founderAllele = Integer.valueOf(words[r+4]);
                    } catch (Exception ex) {}
                    if (founderAllele<0 || founderAllele>=popdata.founderAlleleCount ) {
                        throw new Exception(fName+".hsa, line "+(line+1)+
                                ": invalid founder allele '"+words[r+4]+"'");
                    }
                    if (founderAllele > maxfounderallele) maxfounderallele = founderAllele;
                    double recpos = Double.NaN;
                    try {
                        recpos = Double.valueOf(hsbwords[r]);
                    } catch (Exception ex) {}
                    if (Double.isNaN(recpos)) {
                        throw new Exception(fName+".hsa/hsb, line "+(line+1)+
                                ": number of non-missing items in both files don't match");
                    }
                    hast.addSegment(recpos/100.0, founderAllele); //cM to Morgan
                }
                r++;
            }
            //System.out.println("line="+line+"haplostruct:"+hast);
            hs[chrnum][homol] = hast;
            line++;
            if (line % (popdata.chromCount()*popdata.ploidy) == 0) {
                popdata.getIndiv(iName).setHaploStruct(hs);
            }
        } while(true);
        //if (line != popdata.indivCount()*popdata.chromCount()*popdata.ploidy) {
        //    throw new Exception(fName+".hsa/hsb: number of lines incorrect");
        //}
        //we now allow new individuals in pedigree that don't occur in haplostruct files
        return(maxfounderallele);
    } //readHaploStructFiles

    public static File getParameterFile(String[] args) {
        String filename;
        if (args.length==0) {
            filename = appName+".par";
        } else {
            filename = args[0];
        }
        try {
            File parFile = new File(filename);
            if (parFile.exists() && parFile.isFile()) {
                return parFile;
            }
            else {
                System.out.println(appName+": parameter file "+filename+" not found.");
                return null;
            }
        }
        catch (Exception ex) {
            System.out.println(appName+": parameter file "+filename+" not found.");
            return null;
        }
    } //getParameterFile


    public static String[] getKeyandValue(String line) {
        String[] words = Tools.readWords(line);
        if (words.length==0 || words[0].charAt(0)==commentChar) {
            return new String[] {"",""};
        }
        String key, value="";
        if (words[0].contains("=")) {
            String[] subwords = words[0].split("=",2);
            key = subwords[0];
            if (subwords.length>1) {
                value = subwords[1];
            }
            else if (words.length>1) {
                value = words[1];
            }
        }
        else {
            key = words[0];
            if (words.length>1) {
                if (words[1].equals("=") && words.length>2) {
                    value = words[2];
                } else {
                    if (words[1].charAt(0)=='=') {
                        value = words[1].substring(1);
                    }
                }
            }
        }
        key = key.toUpperCase();
        return new String[] {key, value};
    } //getKeyandValue



    public static boolean isPositiveInteger(String s) {
        try {
            return Integer.parseInt(s)>0;
        } catch (Exception ex) { return false; }
    } //isPositiveInteger

    public static boolean isInteger(String s) {
        try {
            return Integer.parseInt(s)<=Integer.MAX_VALUE;
        } catch (Exception ex) { return false; }
    } //isInteger

    public static boolean isProbability(String s) {
        try {
            double d = Double.parseDouble(s);
            return d>=0.0 && d<=1.0;
        } catch (Exception ex) { return false; }
    } //isProbability

    public static boolean isValidPloidy(String s) {
        try {
            int p = Integer.parseInt(s);
            //return p==2 || p==4;
            return p>0 && p%2==0;
        } catch (Exception ex) { return false; }
    }

    public static boolean isBoolean(String s) {
        s = s.toUpperCase();
        return s.equals("0") || s.equals("1") ||
                s.equals("FALSE") || s.equals("TRUE") ||
                s.equals("NO") || s.equals("YES");
    } //isValidPloidy

    public static boolean toBoolean(String s) {
        s = s.toUpperCase();
        return s.equals("1") || s.equals("TRUE") || s.equals("YES");
    } //toBoolean

    /**
     * canCreateFile: this will delete an existing file with this name,
     * try to create a new one, and delete it again
     * @param fName
     * @return
     */
    public static boolean canCreateFile(String fName) {
        boolean result = false;
        try {
            File test = new File(fName);
            if (test.exists()) {
                test.delete();
            }
            if (test.createNewFile()) {
                result = test.delete();
            }
        } catch (Exception ex) { }
        return result;
    } //canCreateFile

    public static ArrayList<String[]> generateStandardPopulation(String popType,
            int popSize, String missing) {
        ArrayList<String[]> pedList = new ArrayList<String[]>();
        int numlength = Integer.toString(popSize).length();
        DecimalFormat format = new DecimalFormat("0000000000".substring(0, numlength));
        if (popType.equals(strF1)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"F1_"+format.format(i),"P1","P2"});
            }
        } else
        if (popType.equals(strF2)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            pedList.add(new String[] {"F1","P1","P2"});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"F2_"+format.format(i),"F1","F1"});
            }
        } else
        if (popType.equals(strBC)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            pedList.add(new String[] {"F1","P1","P2"});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"BC_"+format.format(i),"P1","F1"});
            }
        } else
        if (popType.equals(strS1)) {
            pedList.add(new String[] {"P1",missing,missing});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"S1_"+format.format(i),"P1","P1"});
            }
        }
        // else no population created
        return pedList;
    } //generateStandardPopulation

    public static void simulateHaploStructs(PopulationData popdata, int maxfounderallele)
            throws Exception {
        for (Individual ind: popdata.getIndividual()) {
            if (ind.getAllHaploStruct()==null) {
                if (ind.isFounder()) {
                    ind.setHaploStruct(Individual.calcFounderHaplostruct(
                            maxfounderallele, popdata));
                    maxfounderallele += popdata.ploidy;
                } else {
                    ind.setHaploStruct(ind.getParents()[0].mating(ind.getParents()[1]));
                }
            }
        }
    } //simulateHaploStructs

    public static void testMeiosis(PopulationData popdata, String fName, double gammaFactor) {
        boolean printPubTables = false; //works only with a specific set of chromosomes and maps,
            //therefore not a user-selectable parameter. 
            //TODO Set to false in public distribution
        try {
            
            //preparation:
            Calendar cal = new GregorianCalendar();
            SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd_HHmmss");
            cal.getTime();
            fName = fName+"_"+sdf.format(cal.getTime());
            popdata.tools.GAMMAFACTOR = gammaFactor;
            double stdev, tvalue, tprob;
            PrintWriter out = new PrintWriter(
                    new BufferedWriter(new FileWriter(fName+".dat")));
            out.println("PedigreeSim test output\t"+fName);
            out.println();
            out.println("test settings:");
            String[] settings = popdata.popdataSettings();
            for (String setting : settings) {
                out.println(setting);
            }
            out.println();
            
            int gamploidy = popdata.ploidy;
            if (popdata.unreducedGametes == 0) 
                gamploidy = gamploidy/2; // normal haploid gametes
            String[] locGenotypeNames = null;
            int[][] locGenotype = null;
            if (gamploidy==1) {
                locGenotypeNames = new String[] {"(0)","(1)"};
            }
            else if (gamploidy == 2) {
                locGenotypeNames = new String[] { "(0,0)","(0,1)","(0,2)","(0,3)",
                    "(1,1)","(1,2)","(1,3)",  "(2,2)","(2,3)", "(3,3)"};
                locGenotype = new int[][] {
                    {0,1,2,3},
                    {1,4,5,6},
                    {2,5,7,8},
                    {3,6,8,9} 
                };
            }
            // for higher ploidy (gamploidy>2) there are too many locus 
            //genotypes, don't count
            
            //print chromosome data and initialize the overall statistics
            out.println("chromosome\tlength\tstart\tend\tcentromere\tprefPairing\tfracQuadrivalents");
            for (Chromosome chr : popdata.getChromosome()) {
                out.print(chr.getChromName()+"\t"+(chr.getLength()*100)+"\t"+
                    (chr.getHeadPos()*100)+"\t"+(chr.getTailPos()*100)+"\t"+(chr.getCentromerePos()*100));
                int mrkCount = chr.getLocus().size();
                    //dosages>1 only possible with double reduction and/or unreduced gametes
                if (gamploidy <= 2) {
                    chr.locGenotypeCount = new int[locGenotypeNames.length][mrkCount];
                }
                chr.recombCountSel = new int[mrkCount][mrkCount];
                if (popdata.ploidy > 2) { 
                    ((TetraploidChromosome)chr).founderHomSel = new int[mrkCount]; 
                }
                chr.founderAlleleCountCum = new int[popdata.founderAlleleCount][mrkCount];
                chr.founderAlleleDoseCum = new int[popdata.founderAlleleCount][mrkCount][5]; //dosage 0..4
                    //dosages>1 only possible with double reduction and/or unreduced gametes
                chr.mapdistCount = new int[mrkCount][mrkCount];
                chr.recombFrCum = new double[mrkCount][mrkCount];
                chr.recombFrSqCum = new double[mrkCount][mrkCount];
                chr.mapdistCum = new double[mrkCount][mrkCount];
                chr.mapdistSqCum = new double[mrkCount][mrkCount];
                chr.recombPoints = new int[(int)(chr.getLength()*4)];
                chr.founderCount = new int[popdata.ploidy];
                if (gamploidy <= 2) {
                    chr.locGenotypeCount = new int[locGenotypeNames.length][mrkCount];
                }
                if (popdata.ploidy > 2) {
                    TetraploidChromosome chr4 = (TetraploidChromosome)chr;
                    out.print("\t"+chr4.getPrefPairingProb()+"\t"+chr4.getFracQuadrivalents());
                    chr4.founderHomFrCum = new double[mrkCount]; //TODO: same
                    chr4.founderHomFrSqCum = new double[mrkCount]; //TODO: same
                    chr4.bivalentCount = 0;
                    chr4.paralQuadrivalentCount = 0;
                    chr4.crossQuadrivalentCount = 0;
                    chr4.exchangeMidFreq = new int[TetraploidChromosome.freqTableLength+1];
                    chr4.exchangeLengthFreq = new int[TetraploidChromosome.freqTableLength+1];
                    chr4.noExchangeLimCount = 0;
                    chr4.oneExchangeLimCount = 0;
                    chr4.twoExchangeLimCount = 0;
                    chr4.exchangeMidSum = 0.0;
                    chr4.exchangeMidSS = 0.0;
                    chr4.exchangeLengthSum = 0.0;
                    chr4.exchangeLengthSS = 0.0;
                    chr4.quadChiasmaSum = 0;
                    chr4.quadChiasmaSS = 0;
                    for (int q=0; q<chr4.quadrivalentConfigCount.length; q++)
                        chr4.quadrivalentConfigCount[q] = 0;
                }
                out.println();
            }
            out.flush();
            
            //Loop over iterations:
            long timertime = System.nanoTime();
            for (int iter=0; iter<popdata.testIter; iter++) {
                System.out.println("Iteration "+iter);
                for (Chromosome chr : popdata.getChromosome()) {
                    int mrkCount = chr.getLocus().size();
                    // initialize the per-iteration statistics:
                    chr.founderAlleleCountSel = new int[popdata.founderAlleleCount][mrkCount];
                    chr.founderAlleleDoseSel = new int[popdata.founderAlleleCount][mrkCount][5]; //dosage 0..4
                        //dosages>1 only possible with double reduction and/or unreduced gametes
                    chr.recombCountSel = new int[mrkCount][mrkCount];
                    if (popdata.ploidy > 2) { 
                        ((TetraploidChromosome)chr).founderHomSel = new int[mrkCount]; 
                    }
                }

                //simulation of gametes and collection of statistics:
                if (popdata.testPrintGametes != 0) {
                    out.println("mei\tgam\tchr\thom\tHaploStruct");
                }
                for (int m=0; m<popdata.testMeioseCount; m++) {
                    ArrayList<Gamete> gam = popdata.getIndiv(0).doMeiosis();
                    int selgam = 0; //gamete 0 is a completely random gamete
                    //int selgam = popdata.tools.rand.nextInt(gam.size()); //alternative
                    if (popdata.testPrintGametes != 0) {
                        //print all gametes or only selected gamete:
                        for (int g=0; g<gam.size(); g++) {
                            if (g==selgam || popdata.testPrintGametes!=1) {
                                for (int c=0; c<popdata.chromCount(); c++) {
                                    for (int hom=0; hom<gamploidy; hom++) {
                                        out.println(m+"\t"+g+"\t"+c+"\t"+hom+"\t"+
                                        gam.get(g).getHaploStruct(c, hom).toString()); 
                                    }    
                                }
                            }    
                        } 
                    }    
                    int g = selgam; // the selected random gamete
                    for (int c=0; c<popdata.chromCount(); c++) {
                        Chromosome chr = popdata.getChrom(c);
                        int locCount = chr.getLocus().size();
                        //count occurrences of locus genotypes:
                        HaploStruct hs0 = gam.get(g).getHaploStruct(c, 0);
                        for (int loc=0; loc<locCount; loc++) {
                            double locuspos = chr.getLocus().get(loc).getPosition();
                            int founderAll0 = hs0.getFounderAt(locuspos);
                            if (gamploidy==1) {
                                chr.locGenotypeCount[founderAll0][loc]++; 
                            }
                            else if (gamploidy==2) {
                                HaploStruct hs1 = gam.get(g).getHaploStruct(c, 1);
                                int founderAll1 = hs1.getFounderAt(locuspos);
                                chr.locGenotypeCount[locGenotype[founderAll0][founderAll1]][loc]++;
                            }  
                            //with higher ploidy there are too many locus genotypes, don't count
                        }

                        //count occurrence of founder alleles and recombinants per haplostruct:
                        for (int p=0; p<gamploidy; p++) {
                            HaploStruct hs = gam.get(g).getHaploStruct(c, p);
                            chr.founderCount[hs.founderCount()-1] ++;
                            int i = hs.segmentCount()-1;
                            if (i>=chr.recombPoints.length) {
                                int[] temp = new int[chr.recombPoints.length];
                                System.arraycopy(chr.recombPoints, 0, temp, 0, temp.length);
                                chr.recombPoints = new int[i+4]; //some extra space
                                System.arraycopy(temp, 0, chr.recombPoints, 0, temp.length);
                            }
                            chr.recombPoints[i] ++;
                            for (int loc=0; loc<locCount; loc++) {
                                double locuspos = chr.getLocus().get(loc).getPosition();
                                int locFounderAllele = hs.getFounderAt(locuspos);
                                chr.founderAlleleCountSel[locFounderAllele][loc]++;
                                for (int loc2=0; loc2<loc; loc2++) {
                                    int loc2FounderAllele =
                                        hs.getFounderAt(chr.getLocus().get(loc2).getPosition());
                                    if (loc2FounderAllele != locFounderAllele) {
                                        chr.recombCountSel[loc][loc2]++;
                                    }
                                }
                            }
                        } //for p

                        //count occurrences of homozygosity for founderalleles:
                        if (gamploidy >= 2) {
                            boolean[] homoz = gam.get(g).homozygous(c,true);
                            for (int loc=0; loc<locCount; loc++) {
                                ((TetraploidChromosome)(chr)).founderHomSel[loc]
                                    += homoz[loc] ? 1 : 0;
                            }
                        }

                        //count the dosages of the founder alleles:
                        int[][] founderdose = gam.get(g).getDosages(c);
                        for (int loc=0; loc<locCount; loc++) {
                            for (int f=0; f<popdata.founderAlleleCount; f++) {
                                chr.founderAlleleDoseSel[f][loc][founderdose[loc][f]] +=1;
                            } 
                        }
                    } // for c
                } // for m

                //output of cumulative results for current iteration:
                if (popdata.testPrintEachIter) {
                    out.println();
                    out.println("Iteration "+iter);
                    out.println();
                    out.println("Cumulative results for ONE RANDOM gamete per meiosis:");
                    int totalhs = gamploidy * popdata.testMeioseCount;

                    for (Chromosome chr: popdata.getChromosome()) {
                        out.println();
                        out.println("Chromosome "+ chr.getChromName());

                        chr.testPrintAlleleCounts(out, chr.founderAlleleCountSel, 
                                totalhs);

                        chr.testPrintAlleleDosageCounts(out, 
                                chr.founderAlleleDoseSel);

                        out.println();
                        out.println("Recombination frequencies:");
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            out.print(chr.getLocus().get(loc).getLocusName());
                            for (int loc2=0; loc2<loc; loc2++)
                                out.print("\t"+(1.0*chr.recombCountSel[loc][loc2]/totalhs));
                            out.println();
                        }
                        chr.printMapHorizontal(out, false);
                        
                        if (popdata.testPrintMapdistances) {
                            out.println();
                            out.println("Inferred map distances:");
                            double recomb, dist;
                            for (int loc=0; loc<chr.getLocus().size(); loc++) {
                                out.print(chr.getLocus().get(loc).getLocusName());
                                for (int loc2=0; loc2<loc; loc2++) {
                                    recomb = 1.0*chr.recombCountSel[loc][loc2]/totalhs;
                                    dist = Tools.recombToMapdist(recomb,popdata.chiasmaInterference);
                                    out.print("\t"+dist);
                                }
                                out.println();
                            }
                            chr.printMapHorizontal(out, false);
                        }

                        if (popdata.ploidy>2) {
                            int total = popdata.testMeioseCount;
                            out.println();
                            out.println("Homozygosity frequencies for founder alleles:");
                            chr.printMapHorizontal(out, true);
                            out.print("freq");
                            for (int loc=0; loc<chr.getLocus().size(); loc++)
                                out.print("\t"+(1.0*((TetraploidChromosome)(chr)).founderHomSel[loc])/total);
                            out.println();
                        }
                    } //for chr
                } //output per iter

                //add the results for this iter (for the one selected gamete per meiosis) to the totals:
                double totalhs = gamploidy * popdata.testMeioseCount;
                for (Chromosome chr: popdata.getChromosome()) {
                    //add the recomb and map distance results for this iter to the totals
                    double recomb, dist;
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        for (int loc2=0; loc2<loc; loc2++) {
                            recomb = 1.0*chr.recombCountSel[loc][loc2]/totalhs;
                            chr.recombFrCum[loc][loc2] += recomb;
                            chr.recombFrSqCum[loc][loc2] += recomb*recomb;
                            dist = Tools.recombToMapdist(recomb,popdata.chiasmaInterference); //not correct for quadrivalents!
                            if (!Double.isNaN(dist)) {
                                chr.mapdistCount[loc][loc2] ++;
                                chr.mapdistCum[loc][loc2] += dist;
                                chr.mapdistSqCum[loc][loc2] += dist*dist;
                            }
                        }
                    }
                    //add the homozygosity frequency for this iter to the totals
                    if (popdata.ploidy>2) {
                        for (int i=0; i<((TetraploidChromosome)(chr)).founderHomSel.length; i++) {
                            double homozFreq =
                                1.0*((TetraploidChromosome)(chr)).founderHomSel[i]
                                    / popdata.testMeioseCount;
                            ((TetraploidChromosome)(chr)).founderHomFrCum[i] +=
                                      homozFreq;
                            ((TetraploidChromosome)(chr)).founderHomFrSqCum[i] +=
                                      homozFreq*homozFreq;
                        }
                    }
                    //add the founder allele counts and dosages for this iter to the totals
                    for (int f=0; f<popdata.founderAlleleCount; f++) {
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            chr.founderAlleleCountCum[f][loc] +=
                                    chr.founderAlleleCountSel[f][loc];
                            for (int d=0; d<5; d++) {
                                chr.founderAlleleDoseCum[f][loc][d] +=
                                        chr.founderAlleleDoseSel[f][loc][d];
                            }
                        }
                    }
                } //for chr
            } //for iter
            timertime = System.nanoTime() - timertime;

            //output accumulated recombination statistics over all iterations:
            int recombfunction=0; //calculated per chromosome; for publication tables we assume
            TDistribution tDistribution = null;
            if (popdata.testIter>1) {
                tDistribution = new TDistribution(popdata.testIter-1);
            }
            out.println();
            out.println("Accumulated statistics over all "+popdata.testIter+" iterations:");
            for (Chromosome chr: popdata.getChromosome()) {
                out.println();
                out.println("Chromosome "+chr.getChromName());

                out.println();
                out.println("Frequency distribution of the number of recombination points");
                int mi = chr.recombPoints.length-1;
                while (mi>=0 && chr.recombPoints[mi]==0) mi--;
                int sum=0, count=0;
                out.println("points\tfrequency");
                for (int i=0; i<=mi; i++) {
                    out.println(i+"\t"+chr.recombPoints[i]);
                    count += chr.recombPoints[i];
                    sum += i*chr.recombPoints[i];
                }
                out.println("count & mean\t"+count+"\t"+(((double)sum)/count));

                out.println();
                out.println("Frequency distribution of the number of different founder alleles");
                out.println("alleles\tfrequency");
                sum=0; count=0;
                for (int i=0; i<popdata.founderAlleleCount; i++) {
                    out.println((i+1)+"\t"+chr.founderCount[i]);
                    count += chr.founderCount[i];
                    sum += (i+1)*chr.founderCount[i];
                }
                out.println("count & mean\t"+count+"\t"+(((double)sum)/count));
                
                double totalhs = popdata.testIter * gamploidy * popdata.testMeioseCount; 
                chr.testPrintAlleleCounts(out, chr.founderAlleleCountCum, 
                        totalhs);
                chr.testPrintAlleleDosageCounts(out, 
                        chr.founderAlleleDoseCum);
                
                if (gamploidy <= 2) {
                    //else too many locus genotypes, not counted
                    out.println();
                    out.println("Frequency distribution of locus genotypes");
                    chr.printMapHorizontal(out, true);
                    for (int i=0; i<locGenotypeNames.length; i++) {
                        out.print(locGenotypeNames[i]);
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            out.print("\t"+chr.locGenotypeCount[i][loc]);
                        }
                        out.println();
                    }
                }    

                int mainconfig = 0; //0: bivalents
                int mainconfigCount = popdata.testIter*popdata.testMeioseCount; //for diploids
                if (popdata.ploidy > 2) {
                    TetraploidChromosome chr4 = (TetraploidChromosome)chr;
                    out.println();
                    out.println("Homozygosity frequencies for founder alleles:");
                    chr4.printMapHorizontal(out, true);
                    out.print("freq");
                    for (int loc=0; loc<chr4.getLocus().size(); loc++)
                        out.print("\t"+(chr4.founderHomFrCum[loc])
                                / popdata.testIter);
                    out.println();
                    out.print("stdev");
                    for (int loc=0; loc<chr4.getLocus().size(); loc++)
                        out.print("\t"+Tools.sampleStDev(
                                chr4.founderHomFrCum[loc], 
                                chr4.founderHomFrSqCum[loc], 
                                popdata.testIter));
                    out.println();

                    out.println();
                    out.println("Total number of bivalents\t" +
                            chr4.bivalentCount);
                    out.println("Total number of parallel quadrivalents\t" +
                            chr4.paralQuadrivalentCount);
                    out.println("Total number of cross-type quadrivalents\t" +
                            chr4.crossQuadrivalentCount);
                    out.println("quadrivalents\tmeiosecount");
                    for (int q=0; q<chr4.quadrivalentConfigCount.length; q++) {
                        out.println(q+"\t"+chr4.quadrivalentConfigCount[q]);
                    }
                    mainconfig = 0; //0: 2 bivalents
                    mainconfigCount = chr4.bivalentCount/2;
                    if (chr4.paralQuadrivalentCount > mainconfigCount) {
                        mainconfig = 1; //1: parallel quadrivalents
                        mainconfigCount = chr4.paralQuadrivalentCount;
                    }
                    if (chr4.crossQuadrivalentCount > mainconfigCount) {
                        mainconfig = 2; //2: cross-type quadrivalents
                        mainconfigCount = chr4.crossQuadrivalentCount;
                    } 
                    
                    out.println();
                    out.println("Data on the chromosome exchange interval of cross-type quadrivalents:");
                    out.println("Number of cross-type quadrivalents where the chromosome exchange interval is:");
                    out.println("not defined (no chiasmata):\t"+
                            chr4.noExchangeLimCount);
                    out.println("defined on one side only:\t" +
                            chr4.oneExchangeLimCount);
                    out.println("defined on both sides:\t" +
                            chr4.twoExchangeLimCount);

                    out.println();
                    out.println("Distribution of midpoints of cross-type quadrivalent exchange intervals (only if chiasmata on both sides)");
                    out.println("interval (M)\t\t\tfrequency");
                    for (int i=0; i<chr4.exchangeMidFreq.length; i++) {
                        out.println((chr4.getHeadPos()+i*(chr4.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t<=X<\t" +
                                (chr4.getHeadPos()+(i+1)*(chr4.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t" + (((double)chr4.exchangeMidFreq[i])/
                                 chr4.twoExchangeLimCount));
                    }
                    double mean =  chr4.exchangeMidSum /
                            chr4.twoExchangeLimCount;
                    out.println("mean\t" + mean);
                    out.println("st.dev.\t" + Tools.sampleStDev(
                            chr4.exchangeMidSum,
                            chr4.exchangeMidSS,
                            chr4.twoExchangeLimCount));
                    
                    out.println();
                    out.println("Distribution of lengths of cross-type quadrivalent exchange intervals (only if chiasmata on both sides)");
                    out.println("interval\t\t\tfrequency");
                    //Note: with popdata.quadriMethod==2 and no chiasma interference the length is always 0: 
                    //the last failed chiasma extends the exchangeLim to the opposing chiasma or chromosome end
                    for (int i=0; i<chr4.exchangeLengthFreq.length; i++) {
                        out.println((chr4.getHeadPos()+i*(chr4.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t<=X<\t" +
                                (chr4.getHeadPos()+(i+1)*(chr4.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t" + (((double)(chr4.exchangeLengthFreq[i]))/
                                 chr4.twoExchangeLimCount));
                    }
                    mean =  chr4.exchangeLengthSum /
                            chr4.twoExchangeLimCount;
                    out.println("mean\t" + mean);
                    out.println("st.dev.\t" + Tools.sampleStDev(
                            chr4.exchangeLengthSum,
                            chr4.exchangeLengthSS,
                            chr4.twoExchangeLimCount));
                } //if ploidy==4

                out.println();
                out.print("Expected recombination frequencies: ");
                recombfunction=0; //bivalents
                if (popdata.ploidy==2 || 
                    mainconfigCount<popdata.testIter*popdata.testMeioseCount ||
                    mainconfig==0) {
                    if (popdata.chiasmaInterference) out.println("Kosambi");
                    else out.println("Haldane");
                }
                else if (mainconfig==1) {
                    recombfunction=1;
                    if (popdata.chiasmaInterference) out.println("Sved_parallel based on Kosambi");
                    else out.println("Sved_parallel based on Haldane");
                }
                else {
                    assert mainconfig==2;
                    recombfunction=2;
                    if (popdata.chiasmaInterference) out.println("Sved_cross-type based on Kosambi");
                    else out.println("Sved_cross-type based on Haldane");
                }
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double dist = chr.getLocus().get(loc).getPosition() -
                                chr.getLocus().get(loc2).getPosition();
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        out.print("\t"+expRec);
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Observed Recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                        out.print("\t"+rec);
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Deviations from expected recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                        double dist = chr.getLocus().get(loc).getPosition() -
                                chr.getLocus().get(loc2).getPosition();
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        out.print("\t"+(rec-expRec));
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Standard Deviations (sample) of recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        out.print("\t"+Tools.sampleStDev(
                                chr.recombFrCum[loc][loc2],
                                chr.recombFrSqCum[loc][loc2],
                                popdata.testIter));
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                if (popdata.testIter>1) {
                    out.println();
                    out.println("Significance of deviations from expected recombination frequencies:");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                            double dist = chr.getLocus().get(loc).getPosition() -
                                    chr.getLocus().get(loc2).getPosition();
                            double expRec = recombfunction==0 ? 
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                    recombfunction==1 ?
                                    Tools.recombBivalentToRecombParallelQuadrivalent(
                                        Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                    Tools.recombBivalentToRecombCrossQuadrivalent(
                                        Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                            stdev= Tools.sampleStDev(
                                    chr.recombFrCum[loc][loc2],
                                    chr.recombFrSqCum[loc][loc2],
                                    popdata.testIter);
                            tvalue = -Math.abs(rec-expRec)/(stdev/Math.sqrt(popdata.testIter));
                            tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                            out.print("\t"+tprob);
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);
                }

                if (popdata.testPrintMapdistances) {
                    out.println();
                    out.println("number of iterations with recombination >= 0.5 :");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            out.print("\t"+(popdata.testIter-chr.mapdistCount[loc][loc2]));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Average inferred map distances (based on bivalent map function):");
                    out.println();
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                            out.print("\t"+meandist);
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Deviations from true distances (based on bivalent map function):");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                            double expDist = chr.getLocus().get(loc).getPosition() -
                                    chr.getLocus().get(loc2).getPosition();
                            out.print("\t"+(meandist-expDist));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Standard Deviations (sample) of map distances:");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            out.print("\t"+Tools.sampleStDev(
                                    chr.mapdistCum[loc][loc2],
                                    chr.mapdistSqCum[loc][loc2],
                                    chr.mapdistCount[loc][loc2]));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    if (popdata.testIter>1) {
                        out.println();
                        out.println("Significance of deviations from true distances (based on bivalent map function):");
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            out.print(chr.getLocus().get(loc).getLocusName());
                            for (int loc2=0; loc2<loc; loc2++) {
                                double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                                stdev= Tools.sampleStDev(
                                        chr.mapdistCum[loc][loc2],
                                        chr.mapdistSqCum[loc][loc2],
                                        chr.mapdistCount[loc][loc2]);
                                double expDist = chr.getLocus().get(loc).getPosition() -
                                        chr.getLocus().get(loc2).getPosition();
                                tvalue = -Math.abs(meandist-expDist)/(stdev/Math.sqrt(chr.mapdistCount[loc][loc2]));
                                if (chr.mapdistCount[loc][loc2]==popdata.testIter)
                                    tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                                else if (chr.mapdistCount[loc][loc2]<2)
                                    tprob = Double.NaN;
                                else tprob = 1.0 - new TDistribution(chr.mapdistCount[loc][loc2]).cumulative(tvalue);
                                out.print("\t"+tprob);
                            }
                            out.println();
                        }
                        chr.printMapHorizontal(out, false);
                    }
                }    

            } // for chr
            out.flush();
            
            if (printPubTables) {
                //Note: this only works with a specific set of
                //chromosomes and maps so printPubTables should generally be set to false
                try {
                    out.println();
                    out.print("Recombination table; expected recombination based on ");
                    if (recombfunction==0) {
                        if (popdata.chiasmaInterference) out.println("Kosambi");
                        else out.println("Haldane");
                    }
                    else if (recombfunction==1) {
                        if (popdata.chiasmaInterference) out.println("Sved_parallel based on Kosambi");
                        else out.println("Sved_parallel based on Haldane");
                    }
                    else {
                        assert recombfunction==2;
                        if (popdata.chiasmaInterference) out.println("Sved_cross-type based on Kosambi");
                        else out.println("Sved_cross-type based on Haldane");
                    }    
                    int[][][][] recloc = new int[][][][] 
                    {
                        { //chrom A
                            { {0,1},{0,2},{0,3},{0,4},{0,6} },
                            { {9,11},{8,12},{7,13},{6,14},{5,15},{0,20} },
                            { {19,20},{18,20},{17,20},{16,20},{14,20} } 
                        }, 
                        { //chrom B
                            { {0,1},{0,2},{0,3},{0,4},{0,5},{0,12} },
                            { {11,13},{10,14},{9,15},{8,16},{7,17},{6,18},{0,24} },
                            { {23,24},{22,24},{21,24},{20,24},{19,24},{12,24} }
                        },
                        { //chrom C
                            { {0,1},{0,2},{0,3},{0,4},{0,5},{0,7},{0,15} },
                            { {14,16},{13,17},{12,18},{11,19},{10,20},{8,22},{7,23},{0,30} },
                            { {29,30},{28,30},{27,30},{26,30},{25,30},{23,30},{15,30} }
                        },
                        { //chrom D
                            { {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,8},{0,17},{0,34} },
                            { {16,18},{15,19},{14,20},{13,21},{12,22},{11,23},{9,25},{8,26},{0,34} },
                            { {33,34},{32,34},{31,34},{30,34},{29,34},{28,34},{26,34},{17,34} }
                        }
                    };
                    out.println("\t\tA\t\tA\t\tA\t\tB\t\tB\t\tB\t\tC\t\tC\t\tC\t\tD\t\tD\t\tD");
                    out.println("interval\trec\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom");
                    for (int line=0; line<9; line++) {
                        double dist = popdata.getChrom(3).getLocus().get(recloc[3][0][line][1]).position;
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        //double expRec = Tools.mapdistToRecomb(dist,popdata.chiasmaInterference);
                        out.print(dist+"\t"+expRec);
                        for (int c=0; c<4; c++) {
                            Chromosome ch = popdata.getChrom(c);
                            for (int col=0; col<3; col++) {
                                if (recloc[c][col].length>line) {
                                    int locA = recloc[c][col][line][0];
                                    int locB = recloc[c][col][line][1];
                                    double rec = ch.recombFrCum[locB][locA]/popdata.testIter;
                                    stdev = Tools.sampleStDev(
                                            ch.recombFrCum[locB][locA],
                                            ch.recombFrSqCum[locB][locA],
                                            popdata.testIter);
                                    tvalue = -Math.abs(rec-expRec)/(stdev/Math.sqrt(popdata.testIter));
                                    tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                                    out.print("\t"+rec+"\t"+tprob);
                                }
                                else out.print("\t\t");
                            }    
                        }
                        out.println();
                    } //for line
                } catch(Exception ex) {
                    System.out.println("Exception in printPubTables: "+ex.getMessage());
                }    
            } //printPubTables
            
            out.println("\nElapsed time\t"+(1.0*timertime/1e9)+"\tseconds");

            out.flush();
        } catch (Exception ex) {
            System.out.println("testMeiosis exception: "+ex.getMessage());
        }
    } //testMeiosis

}