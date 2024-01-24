#####################################################################################################################
# This scripts will be used to generate the embedding vector for a given protein sequence. 
# THe models used are the ones from the Rostlab: https://huggingface.co/Rostlab
# The models are: prot_bert, and prot_t5
# More models will be added in the future. 
# Meng, 01/20/2024 
#####################################################################################################################

from transformers import T5Tokenizer, T5EncoderModel
from transformers import BertModel, BertTokenizer
import torch
import re
import sentencepiece

################################# set up device ######################################
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu' )
print("using device: {}".format(device))


################################ read the protein sequence ############################




################################# prot_bert ##########################################  
tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False )
model = BertModel.from_pretrained("Rostlab/prot_bert")
sequence_Example = "A E T C Z A O"
sequence_Example = re.sub(r"[UZOB]", "X", sequence_Example)
encoded_input = tokenizer(sequence_Example, return_tensors='pt')
embedding_output = model(**encoded_input)
print(embedding_output)

'''
################################# prot_t5 ##########################################
transformer_link = "Rostlab/prot_t5_xl_half_uniref50-enc"
print("loading: {}".format(transformer_link))

model = T5EncoderModel.from_pretrained(transformer_link)
model.full() if device=='cpu' else model.half()
model = model.to(device)
tokenizer = T5Tokenizer.from_pretrained(transformer_link, do_lower_case=False)

sequence_examples = ["proteininfo", "sequence"]

# this will replace all rare or ambiguous amino acids by X and introduce white space between all amino acids

sequence_examples = [" ".join(list(re.sub("[UZOB]", "X", sequence)))  for sequence in sequence_examples ]

# tokenize sequences and pad up to the longest sequence in the batch
ids = tokenizer.batch_encode_plus(sequence_examples, add_special_tokens=True, padding="longest")
input_ids = torch.tensor(ids['input_ids']).to(device)
attention_mask = torch.tensor(ids['attention_mask']).to(device)


# generate embeddings
with torch.no_grad():
  embedding_repr = model(input_ids=input_ids. attention_mask=attention_mask)
  '''