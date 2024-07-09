import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa
import sys
from datetime import datetime
#%pip install pyarrow

#set up data
name=sys.argv[1]
timestamp=datetime.today().strftime('%Y_%m_%d')
data=pd.read_csv(name, delimiter="\t",header=None)
data_str=data.astype('string')
data_pa_table=pa.Table.from_pandas(data)
pq.write_table(data_pa_table,name+'.parquet', version='1.0')
#data.to_parquet('update_isolate_crosswalk.parquet', index=False, engine="pyarrow")
