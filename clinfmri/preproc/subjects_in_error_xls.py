import os
import xlwt


timepoint = "BL"

if timepoint == "BL":
    source_dir = "/neurospin/imagen/BL/processed/spmpreproc_V2"
elif timepoint == "FU2":
    source_dir = "/neurospin/imagen/FU2/processed.new/spmpreproc_V2"
else:
    raise Exception("Timepoint not recognized")
output_xls = "/home/rc244162/{}_spm_preproc_qc.xls".format(timepoint)

subjects_comments = {
    "BL":
        {"000051289725": "image surechantillonnee en z",
         "000074201575": "image accumulee en z",
         "000054713646": "36 coupes en z au lieu de 40",
         "000058908172": "36 coupes en z au lieu de 40",
         "000019767938": "36 coupes en z au lieu de 40",
         "000070675464": "36 coupes en z au lieu de 40",
         "000030775893": "36 coupes en z au lieu de 40",
         "000066715400": "36 coupes en z au lieu de 40"},
    "FU2":
        {"000020401625" : "36 coupes en z au lieu de 40",
         "000028359931" : "36 coupes en z au lieu de 40",
         "000036043675" : "accumulation en z du nifti",
         "000073628015" : "accumulation en z du nifti",
         "000094506022" : "accumulation en z du nifti"}
}

subjects = set()
for item in os.listdir(source_dir):
    if item.isdigit():
        if len(item) == 12:
            subjects.add(item)
for subject, _ in subjects_comments[timepoint].iteritems():
    subjects.add(subject)
subjects = sorted(list(subjects))
# Create an xls object
w = xlwt.Workbook()
ws = w.add_sheet('sheet1')
# Write headers
ws.write(0, 0, "Subject ID")
ws.write(0, 1, "spm_preproc_qc")
ws.write(0, 2, "comment")

# Start filling
row_idx = 1
for subject in subjects:
    ws.write(row_idx, 0, subject)
    if subject in subjects_comments[timepoint]:
        ws.write(row_idx, 1, 1)
        ws.write(row_idx, 2, subjects_comments[timepoint][subject])
    else:
        ws.write(row_idx, 1, 0)
    row_idx += 1

# Save file to disk
w.save(output_xls)
