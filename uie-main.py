"""
@File ：fetch_gene_part.py
@Author ：wu_leanne
@Date ：2022/9/30
@license ：@EndoAngel
@Desc : an example for how to generate the risk level
@Email : wu_leanne@163.com
"""

from paddlenlp import Taskflow
import pandas as pd
from confirm_risk_level.instestines_deal_core import IntestinesObject

uie_model_path = './checkpoint'
schema = ["polyp site", "polyp number", "polyp size", "polyp pathology"]
ie = Taskflow("information_extraction", schema=schema, task_path=uie_model_path)
level_standard = "China"


def main(patient_data_file_path):
    df = pd.read_excel(patient_data_file_path)
    for i in range(len(df)):
        check_id = df.loc[i, "check_id"]
        microscopic_diagnosis = df.loc[i, "microscopic_diagnosis"]
        microscopic_desc = df.loc[i, "microscopic_desc"]
        pathology = df.loc[i, "pathology"]
        entity_list = ie(microscopic_desc)
        intestines_obj = IntestinesObject(check_id, microscopic_desc, microscopic_diagnosis, pathology, level_standard,
                                          entity_list)
        df.loc[i, 'AI-result'] = intestines_obj.level
    df.to_excel(r"./result.xlsx", index=False)
    print(f"Success: Save the result to : ./result.xlsx")


if __name__ == '__main__':
    patient_data = r"./Patient_data.xlsx"
    main(patient_data)
