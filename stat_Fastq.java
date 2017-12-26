/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
/**
 *
 * @author haibo.li
 */
public class stat_Fastq {
        static long totalReads=0;
        static long totalBases=0;
        static long totalQ20=0;
        static long totalQ30=0;
        static long totalGC=0;
        static long totalN=0;
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
	if(args.length<3){
		System.out.println("  Version 1.0.2");
		System.out.println("  Author haibo.li");
		System.out.println("");
		System.out.println("  Usage: java stat_Fastq <33|64> <Output_prefix> <file1.fq|file1.fq.gz> [file2.fq|file2.fq.gz] [...]\n");
		System.exit(0);
	}
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	System.out.println("Begin:\t"+df.format(new Date()));
	int offset=Integer.parseInt(args[0]);
        Read.offset=Integer.parseInt(args[0]);
	String Prefix=(new StringBuffer(args[1])).toString();

	for(int i=2;i<args.length;i++){
		String Reads_fq=args[i];
		File f=new File(Reads_fq);
		if (f.exists()){
			String fileName = f.getName();
			if (fileName.indexOf("gz") == -1){
				String name=fileName.replaceAll("(.fastq|.fq)","");
				stat(f,name,"flat");
				Read.ReadCount=0;
				Read.BaseCount=0;
			}else{
				String name=fileName.replaceAll("(.fastq.gz|.fq.gz)","");
				stat(f,name,"gz");
				Read.ReadCount=0;
				Read.BaseCount=0;
			}

		}else{
			System.out.println(Reads_fq+" file does not exist!");
			System.exit(0);
		}
	}

	String FileName=Prefix+".stat";
	File Basic_Stat =new File(FileName);
	BufferedWriter bw_Stat = new BufferedWriter(new FileWriter(Basic_Stat));
	bw_Stat.write("Sample\tReads\tBases\tAve_length\tQ20(%)\tQ30(%)\tGC(%)\tNppm\n");
	bw_Stat.write(String.format("%s\t%d\t%d\t%.1f\t%.2f\t%.2f\t%.2f\t%d\n",Prefix,totalReads,totalBases,(totalBases*1.0)/totalReads,(totalQ20*100.0)/totalBases,(totalQ30*100.0)/totalBases,(totalGC*100.0)/totalBases,(totalN*1000000)/totalBases));
	bw_Stat.close();
	System.out.println("End:\t"+df.format(new Date()));		
    }
        private static void stat(File f,String name,String format) throws IOException{
        long count[][]=new long[1000][42];
        long actg[][]=new long[1000][20];
        long aveq[]=new long[42];
        int min=1000,max=0;

        File Basic_Info =new File(name+".BasicInfo");
        File GC_content =new File(name+".GCContent");
        File Qual_Dist =new File(name+".QualDist");
        File Read_MeanQual =new File(name+".MeanQual");
        File Base_Qual =new File(name+".qual");
        
        BufferedWriter bw_Basic = new BufferedWriter(new FileWriter(Basic_Info));
        BufferedWriter bw_GC = new BufferedWriter(new FileWriter(GC_content));
        bw_GC.write("#pos\tA\tC\tG\tT\tN\tGCContent\n");
        BufferedWriter bw_Dist = new BufferedWriter(new FileWriter(Qual_Dist));
        BufferedWriter bw_MeanQual = new BufferedWriter(new FileWriter(Read_MeanQual));
        bw_MeanQual.write("#MeanQVofRead\tCount\n");
        BufferedWriter bw_qual = new BufferedWriter(new FileWriter(Base_Qual));

        String title_line,base_line,flag_line,qual_line;
        
        try {
			BufferedReader br;
			if(format == "gz"){
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
			}else{
				br = new BufferedReader(new FileReader(f));
			}
            
            while ((title_line=br.readLine())!=null) {
                base_line=br.readLine();
                flag_line=br.readLine();
                qual_line=br.readLine();
                Read tmp_read=new Read(title_line,base_line,qual_line);
                
                max=(tmp_read.len>max)?tmp_read.len:max;
                min=(tmp_read.len<min)?tmp_read.len:min;
                
                for(int i=0;i<tmp_read.seqBs.length;i++){
                    actg[i][tmp_read.seqBs[i]]++;
                }
                
                for(int i=0;i<tmp_read.qualBs.length;i++){
                    count[i][tmp_read.qualBs[i]]++;
                }
                
                aveq[tmp_read.aveQual]++;
            }
            br.close();
            
            long Q30=0,Q20=0,Q10=0;
            for (int i = 0; i < max; i++) {
                for (int j = 10; j < count[i].length; j++) { Q10+=count[i][j]; }
                for (int j = 20; j < count[i].length; j++) { Q20+=count[i][j]; }
                for (int j = 30; j < count[i].length; j++) { Q30+=count[i][j]; }
            }

            long GC_base=0;
            for (int i = 0; i < max; i++) {
                GC_base+=actg[i][2];
		GC_base+=actg[i][6];
            }

            long N_base=0;
            for (int i = 0; i < max; i++) {
                N_base+=actg[i][13];
            }

		totalReads+=Read.ReadCount;
		totalBases+=Read.BaseCount;
		totalQ20+=Q20;
		totalQ30+=Q30;
		totalGC+=GC_base;
		totalN+=N_base;

            bw_Basic.write(String.format("Total Reads  Count (#):\t%d\n",Read.ReadCount));
            bw_Basic.write(String.format("Total Bases  Count (#):\t%d\n",Read.BaseCount));
            bw_Basic.write(String.format("Average Read Length(bp):\t%.1f\n",Read.BaseCount*1.0/Read.ReadCount));
            bw_Basic.write(String.format("Maxium  Read Length(bp):\t%d\n",max));
            bw_Basic.write(String.format("Minium  Read Length(bp):\t%d\n",min));
            bw_Basic.write(String.format("Q30 Bases Count(bp):\t%d\n",Q30));
            bw_Basic.write(String.format("Q30 Bases Ratio(%%) :\t%.2f\n",Q30*100.0/Read.BaseCount));
            bw_Basic.write(String.format("Q20 Bases Count(bp):\t%d\n",Q20));
            bw_Basic.write(String.format("Q20 Bases Ratio(%%) :\t%.2f\n",Q20*100.0/Read.BaseCount));
            bw_Basic.write(String.format("Q10 Bases Count(bp):\t%d\n",Q10));
            bw_Basic.write(String.format("Q10 Bases Ratio(%%) :\t%.2f\n",Q10*100.0/Read.BaseCount));
            bw_Basic.write(String.format("GC  Bases Count(bp):\t%d\n",GC_base));
            bw_Basic.write(String.format("GC  Bases Ratio(%%) :\t%.2f\n",GC_base*100.0/Read.BaseCount));
            bw_Basic.write(String.format("N   Bases Count(bp):\t%d\n",N_base));
            bw_Basic.write(String.format("N   Bases Ratio(%%) :\t%.2f\n",N_base*1.0/Read.BaseCount));
            
            for (int i = 0; i < max; i++) {
                bw_GC.write((i+1)+"\t"+actg[i][0]+"\t"+actg[i][2]+"\t"+actg[i][6]+"\t"+actg[i][19]+"\t"+actg[i][13]+"\t"+(actg[i][2]+actg[i][6])+"\n");
            }
            
            for (int i = 2; i <= 41; i++) {
                bw_MeanQual.write(i+"\t"+aveq[i]+"\n");
            }
            
            for (int i = 0; i < max; i++) {
                bw_qual.write((i+1)+"");
                for (int j = 2; j <= 41; j++) {
                    bw_qual.write("\t"+count[i][j]);
                }
                bw_qual.write("\n");
            }
            
            int tmp_qual[];
            for (int i = 0; i < max; i++) {
                int read_count=0;
                for (int j = 2; j <= 41; j++) {
                    read_count+=count[i][j];
                }
                tmp_qual=new int[read_count];
                int index=0;
                for (int j = 2; j <= 41; j++) {
                    for (int k = 0; k < count[i][j]; k++) {
                        tmp_qual[index]=j;
                        index++;
                    }
                }
                bw_Dist.write(statQuartile(tmp_qual));
            }
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(stat_Fastq.class.getName()).log(Level.SEVERE, null, ex);
        }
        bw_Basic.close();
        bw_GC.close();
        bw_Dist.close();
        bw_MeanQual.close();
        bw_qual.close();
                
    }
    private static String statQuartile(int arr[]){
        int size=arr.length;
        int mid, mid1, mid3;
        double median, median1, median3, max, min,IQR,tmp_max,tmp_min;
        mid = size / 2;
        median=(size%2 == 0)?(arr[mid-1] + arr[mid])/2.0:arr[mid];
        mid1 = mid / 2;
        mid3 = size%2 == 0  ? (mid + mid1) : (mid+mid1+1);
        median1 = (mid % 2 == 0)?(arr[mid1] + arr[mid1-1])/2.0:arr[mid1];
        median3 = (mid % 2 == 0)?(arr[mid3] + arr[mid3-1])/2.0:arr[mid3];
        IQR = median3 - median1;
	tmp_max = median3 + 1.5*IQR;
	tmp_min = median1 - 1.5*IQR;
 	min = tmp_min < arr[0]?arr[0]:tmp_min;
	max = tmp_max > arr[size-1]?arr[size-1]:tmp_max;
        return(String.format("%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",min,median1,median,median3,max));
    } 
    
}

class Read{
    static int offset=33;
    static long ReadCount=0;
    static long BaseCount=0;
	static int len=0,minLen=10000,maxLen=0;
    static int Q30=0,Q20=0,Q10=0;
    static int GC=0,N=0;
    String ID,seq,qual;
	String name="+";
    int aveQual=0;
    int qualBs[];
	int seqBs[];
	
    public Read(String nameString,String seqString,String qualString){
        String title_array[]=nameString.split("\\s");
        ID=title_array[0];
        seq=seqString; 
        qual=qualString;
        len=seq.length();
        ReadCount++;
        BaseCount+=len;
        getSeqArray();
        getQualArray();
        aveQual=getAveQual();
        char ch[]=seq.toCharArray();
        if(minLen>len){
            minLen=len;
        }
        if(maxLen<len){
            maxLen=len;
        }
        ReadCount++;
        BaseCount+=len;
        for(int i=0;i<ch.length;i++){
            if(ch[i]=='G'||ch[i]=='C'){
                GC+=1;
            }
        }
        for(int i=0;i<ch.length;i++){
            if(ch[i]=='N'){
                N+=1;
            }
        }
        Q30+=getQual(30);
        Q20+=getQual(20);
        Q10+=getQual(10);
	}
    void getSeqArray(){
        byte tmp_byte[]= seq.getBytes();
        seqBs=new int[tmp_byte.length];
        for (int i = 0; i < tmp_byte.length; i++) {
           seqBs[i]=(int) tmp_byte[i]-65;
        }        
    }
    void getQualArray(){
        byte tmp_byte[]= qual.getBytes();
        qualBs=new int[tmp_byte.length];
        for (int i = 0; i < tmp_byte.length; i++) {
            int tmp_qual=(int) tmp_byte[i]-offset;  
            qualBs[i]=(tmp_qual>=41)?41:tmp_qual;
            //qualBs[i]=(int) tmp_byte[i]-offset;
        }
    }
    int getAveQual(){
        int sumQual=0;
        for (int i = 0; i < qualBs.length; i++) {
            sumQual+=qualBs[i];
        }
        return aveQual = sumQual/qualBs.length;
    }
	    int getQual(int threshold){
        int count=0;
        byte tmp_byte[]= qual.getBytes();
        for (int i = 0; i < tmp_byte.length; i++) {
            int tmp_qual=(int) tmp_byte[i]-offset;  
            if (tmp_qual>=threshold) {
                count+=1;
            }
        }
        return count;
    }
}
