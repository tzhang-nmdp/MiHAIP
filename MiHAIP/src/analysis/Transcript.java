package analysis;

import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by wwang on 9/5/16.
 */
public class Transcript {
    public int chrome;
    public int pos;
    String geneID;
    String ref;
    String alt;
    String transcriptID;
    Patient patient;
    //true if the transcript is on the minus strand (CDS is the reverse complement of the VCF bases)
    boolean minusStrand = false;
    public HashMap<Integer, String> change = new HashMap<>();

    public Transcript(String transcriptID){
        this.transcriptID = transcriptID;
    }

    public void setInfo(Transcript t){
        this.chrome = t.chrome;
        this.pos = t.pos;
        this.ref = t.ref;
        this.alt = t.alt;
        this.patient = t.patient;
        this.minusStrand = t.minusStrand;
    }

    /**
     * Derive the strand from the SnpEff HGVS.c field (e.g. c.380G>C) and the VCF ref base.
     * HGVS.c ref base == VCF ref -> plus strand; == complement of VCF ref -> minus strand.
     * Anything that is not a single base substitution keeps the current strand.
     */
    public void setStrandFromHgvsC(String hgvsC){
        if(hgvsC == null || ref == null || ref.length() != 1){
            return;
        }
        Matcher m = Pattern.compile("c\\.\\d+([ACGT])>[ACGT]").matcher(hgvsC);
        if(!m.find()){
            return;
        }
        String cRef = m.group(1);
        String vcfRef = ref.toUpperCase();
        if(cRef.equals(vcfRef)){
            minusStrand = false;
        }else if(cRef.equals(complement(vcfRef))){
            minusStrand = true;
        }
    }

    public static String complement(String s){
        StringBuilder sb = new StringBuilder();
        for(char c : s.toCharArray()){
            switch (c){
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                case 'a': sb.append('t'); break;
                case 't': sb.append('a'); break;
                case 'c': sb.append('g'); break;
                case 'g': sb.append('c'); break;
                default: sb.append(c);
            }
        }
        return sb.toString();
    }

    public void addPos(int pos, String NT){
        change.put(pos, NT);
    }

    public void setGeneID(String geneID){
        this.geneID = geneID;
    }

    public void setNT(String ref, String alt){
        this.ref = ref;
        this.alt = alt;
    }
    public String getAlt(){
        return alt;
    }
    public void setChrome(int i){
        chrome = i;
    }
    public void setChromePos(int i){
        pos = i;
    }
    public void setPatient(Patient p){
    	this.patient = p;
    }
    public String getInfo(){
    	if(chrome == 24){
    		return ">chr"+"X" + "|" + pos + "|"+ geneID + "|" + patient.getID();
    	}else if(chrome == 25){
    		return ">chr"+"Y" + "|" + pos + "|"+ geneID+ "|" + patient.getID();
    	}else if(chrome == 26){
    		return ">chr"+"M" + "|" + pos + "|"+ geneID+ "|" + patient.getID();
    	}else{
    		return ">chr"+chrome + "|" + pos + "|"+ geneID+ "|" + patient.getID();
    	}
        
    }
    public String getChrome(){
    	if(chrome == 24){
    		return "chrX";
    	}else if (chrome == 25){
    		return "chrY";
    	}else{
    		return "chr"+chrome;
    	}
        
    }
}
