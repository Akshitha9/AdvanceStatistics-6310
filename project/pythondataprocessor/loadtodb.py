import sqlite3
import pandas as pd
from tqdm import tqdm
import os

dropandcreate = 'Y'
loadbasetables = 'Y'
dbinfo = 'db/project.db'


def main():
    if dropandcreate == 'Y':
        with sqlite3.connect(dbinfo) as conn:
            cur = conn.cursor()
            #cur.execute("DROP TABLE data_mrna_seq_v2_rsem_zscores_ref_all_samples")
            cur.execute("DROP TABLE patients")
            #cur.execute("delete from genevalues")
            cur.execute("delete from genelist")
            conn.commit()

    if loadbasetables == 'Y':
        df = pd.read_csv('sourcedata/patientsdata.csv')
        df1 = pd.read_csv('sourcedata/data_mrna_seq_v2_rsem_zscores_ref_all_samples.csv')

        with sqlite3.connect(dbinfo) as con:
            df.to_sql("patients", con)
            con.commit()
            cur = con.cursor()
            cur.execute(r'''select PATIENT_ID from patients where VITAL_STATUS is null or age is null or sex is null or ALCOHOL_CONSUMPTION is null or TOBACCO_SMOKING_HISTORY is null
                        union
                        select patient_id  from patients where VITAL_STATUS ='NA' or age ='NA' or sex ='NA' or ALCOHOL_CONSUMPTION ='NA' or TOBACCO_SMOKING_HISTORY ='NA'
                        union
                        select patient_id from  patients where ALCOHOL_CONSUMPTION  IN ('Alcohol consumption history not available','Consumed alcohol in the past, but currently a non-drinker')
                         union
                        select patient_id from patients where TOBACCO_SMOKING_HISTORY  IN ('Current reformed smoker, more than 15 years','Current reformed smoker within past 15 years','Current reformed smoker, years unknown','Smoking history not available')
                        ''')
            
            data = cur.fetchall()
            noteligible = [ i[0] for i in data]
            cur.execute('delete from patients where VITAL_STATUS is null or age is null or sex is null or ALCOHOL_CONSUMPTION is null or TOBACCO_SMOKING_HISTORY is null')
            cur.execute("delete from patients where VITAL_STATUS ='NA' or age ='NA' or sex ='NA' or ALCOHOL_CONSUMPTION ='NA' or TOBACCO_SMOKING_HISTORY ='NA'")
            cur.execute("delete from patients where ALCOHOL_CONSUMPTION  IN ('Alcohol consumption history not available','Consumed alcohol in the past, but currently a non-drinker') ")
            cur.execute("delete from patients where TOBACCO_SMOKING_HISTORY  IN ('Current reformed smoker, more than 15 years','Current reformed smoker within past 15 years','Current reformed smoker, years unknown','Smoking history not available') ")
            cur.close()
            con.commit()
            
            return df1,noteligible
        
def generafiles(df1,noteligible):
    conn = sqlite3.connect(dbinfo)
    df = pd.read_sql_query(r'''select PATIENT_ID,VITAL_STATUS,age,sex,case ALCOHOL_CONSUMPTION when 'Alcohol consumption equal to or less than 2 drinks per day for men and 1 drink or less per day for women' then 'Drinker'
                                                                when 'Alcohol consumption more than 2 drinks per day for men and more than 1 drink per day for women' then 'Drinker' 
                                                                else 'Non-Drinker'
                                                                end drinkingstatus
                                    ,case TOBACCO_SMOKING_HISTORY when 'Current smoker: Includes daily and non-daily smokers' then 'Smoker'
                                     else  'Non-Smoker'
                                    end smokingstatus
                                from patients order by age''',conn)
    df.to_csv(f'finalpatients.csv',index=False)
    
    patientlist = df['PATIENT_ID'].tolist()
    print(patientlist)
    #drop excluded patients
    for i in noteligible:
        df1.drop(i, inplace=True, axis=1)
    
    genecollist = df1.columns.tolist()
    sortlist = genecollist[0:2]+patientlist
    #drop data where gene value is same for all patients
    metadf = df1[genecollist[0:2]]
    patients =  df1[patientlist]
    patients = patients[patients.diff(axis=1).fillna(0).ne(0).any(axis=1)]
    fineldf = pd.concat([metadf, patients], axis=1)
    #drop NA data
    fineldf = fineldf.dropna()
    print(sortlist)
    #df1 = df1[sortlist]
    
    fineldf.to_csv('finalgenevalues.csv')
    

    #df2 = pd.read_sql_query('select * from genevalues LIMIT 5',conn)
    #df2.to_csv(f'finalgenevalues1.csv',index=False)
    
    #pd.read_csv('finalgenevalues1.csv', header=None).T.to_csv('finalgenevalues.csv', header=False, index=False)
    #os.remove('finalgenevalues1.csv')


if __name__ == "__main__":
    df1,noteligible = main()
    generafiles(df1,noteligible)
    
#cd /Users/akshitha/Documents/GitHub/AdvanceStatistics-6310/project/pythondataprocessor 
#source env/bin/activate
#python loadtodb.py
