## Getting started

1. AS system requires Python 3.8+.  You can install requirements `pip`:

   ```
   pip install -r requirements.txt
   ```

2. Download UIE model:

   https://endoangel-info.com/checkpoint.zip

   **When the file is completely downloaded, extract the file to its root directory.**

   - The file directory is as follows

     ```
     ├─checkpoint
     │  │  model_config.json
     │  │  model_state.pdparams
     │  │  special_tokens_map.json
     │  │  tokenizer_config.json
     │  │  vocab.txt
     │  │
     │  └─static
     │          inference.pdiparams
     │          inference.pdiparams.info
     │          inference.pdmodel
     │
     └─confirm_risk_level
     │   │  ai_model.py
     │   │  instestines_deal_core.py
     │   │  pathlogy_model.py
     │   │  risk_level_model.py
     │   │  __init__.py
     │  Patient_data.xlsx
     │  README.md
     │  requirements.txt
     │  uie-main.py
     ```

3. Below is an example for classifying patient risk level and assigning the surveillance interval

   ```
   python uie-main.py
   ```

      > `Patient_data.xlsx`  is  the patient data
      >
      > `result.xlsx` is the risk level result