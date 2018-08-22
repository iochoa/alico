//
//  compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <stdbool.h>
#include <algorithm>
#include "sam_block.h"
#include "read_compression.h"

#include "IO/SAM/SAMRecord.h"
#include "QualCodec/QualEncoder.h"
#include "QualCodec/QualDecoder.h"
#include "Common/Exceptions.h"
#include "IO/CQ/CQFile.h"
#include "IO/File.h"

char idx_convert(int i) {
  if(i==0) {return 'A';}
  else if(i==1) {return 'T';}
  else if(i==2) {return 'C';}
  else if(i==3) {return 'G';}
  return 'N';
}

void reconstruct_ref(void *thread_info) {
  int LCC;
  int context_idx,t;
  uint32_t tag;
  uint64_t previous;
  struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
  Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->frefcom);
  char tt1[3],tt2[3];

  char prefix[2]; prefix[0]='F'; prefix[1]='F';
  uint64_t context[25][6];

  for(int i=0; i<25; i++) {
    for(int j=0; j<6; j++) {
      if(j==5) {context[i][j]=5;}
      else {context[i][j]=1;}
    }
  }

  char head[1024];
  char num[1024];
  char tmp;
  uint32_t total;

  while(!feof(info->fsinchr)) {
    fgets(head, 1024, info->fsinchr);
    fgets(num, 1024, info->fsinchr);
    if(prefix[0]=='F') {
    	LCC=atoi(num);
    	fgets(tt1, 3, info->fsinchr);
    	fgets(tt2, 3, info->fsinchr);
    }
    fscanf(info->fsinchr, "%s\n", num);
    fprintf(info->fref, "%s", head);
    //printf("total num: %s\n", num);
    total = atoi(num);
    for(int i = 0; i<total; i++) {
      if(prefix[0]=='F') {
        fputc(tt1[0],info->fref);
        prefix[0]=tt1[0];
      }
      else if(prefix[1]=='F') {
        fputc(tt2[0],info->fref);
        prefix[1]=tt2[0];
      }
      else {
      //printf("%c %c\n",prefix[0],prefix[1]);
        context_idx = context_index(prefix[0])*5+context_index(prefix[1]);

	//for(int tt=0;tt<6;tt++) {printf("%lu ",context[context_idx][tt]);}
	//printf("\n");

        tag = arithmetic_get_symbol_range(as, ((uint32_t)(context[context_idx][5])));
       //printf("%u\n",tag);
        previous = 0; t=0;
        while(previous<=tag) {
          previous += context[context_idx][t];
          t++;
        }
        t--;
        tmp = idx_convert(t);
        fputc(tmp,info->fref);
        if((i+1)%LCC==0) {fputc('\n',info->fref);}
        arithmetic_decoder_step(as, ((uint32_t)(previous-context[context_idx][t])), ((uint32_t)previous), ((uint32_t)context[context_idx][5]));
        context[context_idx][t]++; context[context_idx][5]++;
	prefix[0]=prefix[1]; prefix[1]=tmp;
      }
    }
    if(total%LCC!=0) {fputc('\n',info->fref);}
  }

  free_os_stream(as->ios);
  free(as);

  return;
}

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs, bool compressing){
    
    char foo[] = "CIGAR";
    
    int32_t i = 0;
    
    switch (print_mode) {
        case 0:
            fprintf(fs, "%s\t", sline->ID);
            fprintf(fs, "%d\t", sline->flag);
            fprintf(fs, "%s\t", sline->rname);
            fprintf(fs, "%d\t", sline->pos);
            fprintf(fs, "%d\t", sline->mapq);
            fprintf(fs, "%s\t", sline->cigar);
            fprintf(fs, "%s\t", sline->rnext); 
            fprintf(fs, "%d\t", sline->pnext);
            fprintf(fs, "%d\t", sline->tlen);
            fprintf(fs, "%s\t", sline->read);
            // need to re-reverse quality scores
            if ((sline->flag & 16) == 16 && !compressing) {
		//printf("sline print line is: %d\n", sline->readLength);
                for (i = sline->readLength - 1; i >= 0; --i)
                    fputc(sline->quals[i], fs);
                fputc('\t', fs);
            } else {
                fprintf(fs, "%s\t", sline->quals);
            }
            fprintf(fs, "%s", sline->aux);  
            fputc('\n', fs); 
            break;
            
        default:
            break;
    }
    return 0;
}

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs) {
    return print_line(sline, print_mode, fs, false);
}


static void build_SAMRecord_for_calq(sam_block samBlock, calq::SAMRecord * const samRecord)
{
    // POS uint32_t = 
    samRecord->pos = samBlock->reads->lines->pos;
    // CIGAR std::string = char*
    samRecord->cigar = samBlock->reads->lines->cigar;
    // SEQ std::string = char*
    // TO-DO: Is seq = read? 
    samRecord->seq = samBlock->reads->lines->read;
    // QUAL std::string = symbol_t*
    for (int i = 0; i <= samBlock->QVs->qv_lines->columns; i++) {
        samRecord->qual += samBlock->QVs->qv_lines->data[i];
    }
}

static void extract_SAM_data_for_calq(sam_block samBlock, uint32_t * const pos, std::string * const cigar, std::string * const seq, std::string * qual) {
    *pos = samBlock->reads->lines->pos;
    *cigar = samBlock->reads->lines->cigar;
    *seq = samBlock->reads->lines->read;
    for (int i = 0; i < samBlock->QVs->qv_lines->columns; i++) {
        *qual += (samBlock->QVs->qv_lines->data[i] + char(33));
    }
    if (samBlock->reads->lines->invFlag & 16) {
        std::reverse(qual->begin(), qual->end());
    }
    printf("%s\n", qual->c_str());
    printf("%s\n", cigar->c_str());
}

int compress_line(Arithmetic_stream as, sam_block samBlock, FILE *funmapped, uint8_t lossiness, int calq, Arithmetic_stream as1, uint64_t context[25][6], char* prefix, FILE * fsinchr, std::vector<calq::SAMRecord> &samRecords)  {
    try {
         /*if (calq) {
            //Max33 Phred+33 [0,93]

                     //}*/
        
        static bool unmapped_reads = false;
        uint8_t chr_change;
        // Load the data from the file
            printf("hi");

        if(load_sam_line(samBlock)){

            return 0;
    }
        // If read is unmapped and reference name is *, we assume that all the remaining
        // lines are unmapped and have reference name *.
        // If the read is unmapped but has a position/reference name, we simply use that
        // information to compress it.
        if ( (samBlock->reads->lines->invFlag & 4) == 4 && *samBlock->rnames->rnames[0] == '*'){
            read_line line = samBlock->reads->lines; 

            fprintf(funmapped, "%s\t", *samBlock->IDs->IDs);
            fprintf(funmapped, "%d\t", line->invFlag);
            fprintf(funmapped, "%s\t", *samBlock->rnames->rnames);
            fprintf(funmapped, "%d\t", line->pos);
            fprintf(funmapped, "%d\t", *samBlock->mapq->mapq);
            fprintf(funmapped, "%s\t", line->cigar);
            fprintf(funmapped, "%s\t", *samBlock->rnext->rnext); 
            fprintf(funmapped, "%d\t", *samBlock->pnext->pnext);
            fprintf(funmapped, "%d\t", *samBlock->tlen->tlen);
            fprintf(funmapped, "%s\t", line->read);
            int32_t i = 0;
            uint32_t j = 0;
            qv_line_t qline = *samBlock->QVs->qv_lines;
            if ((line->invFlag & 16) == 16) {
                for (j = qline.columns; j > 0; j--) {
                    fputc(qline.data[j] + 33, funmapped);
                }
            } else {
                for (j = 0; j < qline.columns; i++) {
                    fputc(qline.data[j] + 33, funmapped);
                }
            }
            fputc('\t', funmapped);
            for (i = 0; i < samBlock->aux->aux_cnt; i++) {
                fprintf(funmapped, "%s", samBlock->aux->aux_str[i]);  
                if (i < samBlock->aux->aux_cnt - 1) fputc('\t', funmapped);
            }
            fputc('\n', funmapped); 
            return 1;
        } else {
            if (unmapped_reads) {
                fprintf(stderr, "compress_line error: There is a mapped read following a read that is not mapped to any chromosome. This probably means that the sam file is not correctly sorted.\n");
                exit(1);
            }
        }
            
        // Compress sam line
        
        chr_change = compress_rname(as, samBlock->rnames->models, *samBlock->rnames->rnames);
    
        if (chr_change == 1){
            // Store Ref sequence in memory
            store_reference_in_memory_com(samBlock->fref, as1, context, prefix, fsinchr);
            // Reset cumsumP
            cumsumP = 0;
            memset(snpInRef, 0, MAX_BP_CHR);
        }
        compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs);

        compress_mapq(as, samBlock->mapq->models, *samBlock->mapq->mapq);

        compress_rnext(as, samBlock->rnext->models, *samBlock->rnext->rnext);

        compress_read(as, samBlock->reads->models, samBlock->reads->lines, chr_change);
        
        compress_cigar(as, samBlock->reads->models, samBlock->reads->lines->cigar, samBlock->reads->lines->cigarFlags);

        compress_tlen(as, samBlock->tlen->models, *samBlock->tlen->tlen);

        compress_pnext_raw(as, samBlock->pnext->models,  samBlock->reads->lines->pos, *samBlock->pnext->pnext);
        
        //idoia
        //compress_aux(as, samBlock->aux->models, samBlock->aux->aux_str, samBlock->aux->aux_cnt, samBlock->aux);
        compress_aux_idoia(as, samBlock->aux->models, samBlock->aux->aux_str, samBlock->aux->aux_cnt, samBlock->aux);
        //idoia
        if (calq) {
            //fprintf(stdout, "line compression w/ calq\n");
            uint32_t pos = 0;
            std::string cigar = "";
            std::string seq = "";
            std::string qual = "";
            extract_SAM_data_for_calq(samBlock, &pos, &cigar, &seq, &qual);
            calq::SAMRecord samRecord(pos, cigar, seq, qual);
            samRecords.push_back(samRecord);
        } else {
            if (lossiness == LOSSY) {
                QVs_compress(as, samBlock->QVs, samBlock->QVs->qArray);
            } else {
                QVs_compress_lossless(as, samBlock->QVs->model, samBlock->QVs->qv_lines);
            }
        }

    }
    catch (const calq::ErrorException& e) {
        std::cerr << "CALQ error";
        if (strlen(e.what()) > 0) {
            std::cerr << ": " << e.what();
        }
        std::cerr << std::endl;
        abort();
    }
    catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        abort();
    }

    return 1;
}

int decompress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness, uint8_t calqmode, std::vector<calq::SAMRecord> &samRecords) {
    
    int32_t chr_change = 0;
    
    uint32_t decompression_flag = 0;
    
    struct sam_line_t sline;
    
    //This is only for fixed length? i think so.
    sline.readLength = samBlock->read_length;
    //sline.readLength = samBlock->reads->models->read_length;
    //printf("sline read length is: %d\n", sline.readLength); 
    //printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
        
    chr_change = decompress_rname(as, samBlock->rnames->models, sline.rname);
    
    if (chr_change == -1)
        return 0;
        
    if (chr_change == 1){
            
        //printf("Chromosome %d decompressed.\n", ++chrCtr);
            
        // Store Ref sequence in memory
        store_reference_in_memory(samBlock->fref);
            
        // reset cumsumP
        cumsumP = 0;

        memset(snpInRef, 0, MAX_BP_CHR);

    }

    decompress_id(as, samBlock->IDs->models, sline.ID);

    decompress_mapq(as, samBlock->mapq->models, &sline.mapq);

    decompress_rnext(as, samBlock->rnext->models, sline.rnext); 

    decompression_flag = decompress_read(as,samBlock, chr_change, &sline);
    
    decompress_cigar(as, samBlock, &sline);

    decompress_tlen(as, samBlock->tlen->models, &sline.tlen);

    decompress_pnext(as, samBlock->pnext->models, sline.pos, sline.tlen, samBlock->read_length, &sline.pnext, sline.rnext[0] != '=', NULL);

    //idoia
    //decompress_aux(as, samBlock->aux, sline.aux);
    decompress_aux_idoia(as, samBlock->aux, sline.aux);
    //idoia
    if (calqmode) {
        //fprintf(stdout, "line compression w/ calq\n");
        uint32_t pos = 0;
        std::string cigar = "";
        std::string seq = "";
        std::string qual = "";
        extract_SAM_data_for_calq(samBlock, &pos, &cigar, &seq, &qual);
        calq::SAMRecord samRecord(pos, cigar, seq, qual);
        samRecords.push_back(samRecord);
        //build_SAMRecord_for_calq(samBlock, &samRecord);
    } else {
   
        if (lossiness == LOSSY) {
            QVs_decompress(as, samBlock->QVs, decompression_flag, sline.quals);
        }
        else
            QVs_decompress_lossless(as, samBlock->QVs, decompression_flag, sline.quals, (int)strlen(sline.read));

    }
    sline.readLength = samBlock->reads->models->read_length;
    // printf("sline read length before printing is: %d\n", sline.readLength); 
    print_line(&sline, 0, samBlock->fs);

    return 1;
}

int compress_most_common_list(Arithmetic_stream as, aux_block aux)
{
    uint8_t n,i,l,k;

    aux_models models = aux->models;
    n = aux->most_common_size;
    compress_uint8t(as,models->most_common_list[0],n);

    for (i=0;i<n;i++) {
        l = strlen(aux->most_common[i]);
        compress_uint8t(as,models->most_common_list[0],l);
        for(k=0;k<l;k++)
            compress_uint8t(as,models->most_common_list[0],aux->most_common[i][k]);
    }

    return 1;
}


int decompress_most_common_list(Arithmetic_stream as, aux_block aux)
{
    uint8_t n,i,l,k;

    aux_models models = aux->models;
    n = decompress_uint8t(as,models->most_common_list[0]);

    aux->most_common_size = n;

    char buffer[256];
    for (i=0;i<n;i++) {
        l = decompress_uint8t(as,models->most_common_list[0]);
        for(k=0;k<l;k++)
            buffer[k]=decompress_uint8t(as,models->most_common_list[0]);
        buffer[k]='\0';

        strcpy(aux->most_common[i],buffer);
        //printf("%d -> %s\n",i,aux->most_common[i]);
    }
    return 1;
}

void free_sam_reads_compress(read_block reads){
    free(reads->lines->cigar);
    free(reads->lines->edits);
    free(reads->lines->read);
    free(reads->lines);
    free(reads->models->flag[0]->counts-1);
    free(reads->models->flag[0]);
    free(reads->models->flag);
    for (int i = 0; i < 4; i++){
    	free(reads->models->pos_alpha[i]->counts-1);
        free(reads->models->pos_alpha[i]);
    }
    free(reads->models->pos_alpha);
    free(reads->models->pos[0]->alphabet);
    free(reads->models->pos[0]->counts);
    free(reads->models->pos[0]->alphaExist);
    free(reads->models->pos[0]->alphaMap);
    free(reads->models->pos[0]);
    free(reads->models->pos);
    for (int i = 0; i < 256; i++){
        free(reads->models->match[i]->counts-1);
        free(reads->models->match[i]);
    }
    free(reads->models->match);
    free(reads->models->snps[0]->counts-1);
    free(reads->models->snps[0]);
    free(reads->models->snps);
    free(reads->models->indels[0]->counts-1);
    free(reads->models->indels[0]);
    free(reads->models->indels);
    for (int i = 0; i < 0xffff; i++){
    	free(reads->models->var[i]->counts-1);
    	free(reads->models->var[i]);
    }
    free(reads->models->var);
    for (int i = 0; i < 6; i++){
        free(reads->models->chars[i]->counts-1);
        free(reads->models->chars[i]);
    }
    free(reads->models->chars);
    free(reads->models->cigar[0]->counts-1);
    free(reads->models->cigar[0]);
    free(reads->models->cigar);
    free(reads->models->cigarFlags[0]->counts-1);
    free(reads->models->cigarFlags[0]);
    free(reads->models->cigarFlags);
    for (int i = 0; i < 4; i++){
        free(reads->models->rlength[i]->counts-1);
        free(reads->models->rlength[i]);
    }
    free(reads->models->rlength);
    free(reads->models);
    free(reads); 
}

void free_sam_ids_compress(id_block IDs){
    free(IDs->IDs[0]);
    free(IDs->IDs);
    for (int i = 0; i < MAX_NUMBER_TOKENS_ID; i++){
        free(IDs->models->alpha_len[i]->counts-1);
        free(IDs->models->alpha_len[i]);
        free(IDs->models->alpha_value[i]->counts-1);
        free(IDs->models->alpha_value[i]);
        free(IDs->models->chars[i]->counts-1);
        free(IDs->models->chars[i]);
        free(IDs->models->delta[i]->counts-1);
        free(IDs->models->delta[i]);
        free(IDs->models->zero_run[i]->counts-1);
        free(IDs->models->zero_run[i]);
        free(IDs->models->token_type[i]->counts-1);
        free(IDs->models->token_type[i]);
        for (int j = 0; j < 4; j++) {
            free(IDs->models->integer[i*4 + j]->counts-1);
            free(IDs->models->integer[i*4 + j]);
        }
    }
    free(IDs->models->alpha_len);
    free(IDs->models->alpha_value);
    free(IDs->models->chars);
    free(IDs->models->integer);
    free(IDs->models->delta);
    free(IDs->models->zero_run);
    free(IDs->models->token_type);
    free(IDs->models);
    free(IDs);
}

void free_sam_block_compress(sam_block samBlock){
    for (int i = 0; i < 256*4; i++){
        free(samBlock->codebook_model[i]->counts-1);
        free(samBlock->codebook_model[i]);
    }
    free(samBlock->codebook_model);
    free_sam_reads_compress(samBlock->reads);
    free_distortion_matrix(samBlock->QVs->dist);
    free_alphabet(samBlock->QVs->alphabet);
    free(samBlock->QVs->qv_lines[0].data);
    free(samBlock->QVs->qv_lines);
    uint32_t model_idx;
    for (uint32_t i = 0; i < samBlock->QVs->columns; i++){
        for (int j = 0; j < 2*(QV_ALPHABET_SIZE + 1); j++){
            model_idx = get_qv_model_index(i, j);
            free(samBlock->QVs->model[model_idx]->counts);
            free(samBlock->QVs->model[model_idx]);
        }
    }
    free(samBlock->QVs->model);
    free_sam_ids_compress(samBlock->IDs);
    for (int i = 0; i < MAX_AUX_FIELDS; i++)
        free(samBlock->aux->aux_str[i]);
    free(samBlock->aux->aux_str);
    for (int i = 0; i < MOST_COMMON_LIST_SIZE; i++)
        free(samBlock->aux->most_common[i]);
    free(samBlock->aux->most_common);
    free(samBlock->aux->models->qAux[0]->counts-1);
    free(samBlock->aux->models->qAux[0]);
    free(samBlock->aux->models->qAux);
    free(samBlock->aux->models->tagtypeLUTflag[0]->counts-1);
    free(samBlock->aux->models->tagtypeLUTflag[0]);
    free(samBlock->aux->models->tagtypeLUTflag);
    free(samBlock->aux->models->typeLUTflag[0]->counts-1);
    free(samBlock->aux->models->typeLUTflag[0]);
    free(samBlock->aux->models->typeLUTflag);
    free(samBlock->aux->models->tagtypeLUT[0]->counts-1);
    free(samBlock->aux->models->tagtypeLUT[0]);
    free(samBlock->aux->models->tagtypeLUT);
    free(samBlock->aux->models->tag[0]->counts-1);
    free(samBlock->aux->models->tag[0]);
    free(samBlock->aux->models->tag[1]->counts-1);
    free(samBlock->aux->models->tag[1]);
    free(samBlock->aux->models->tag);
    free(samBlock->aux->models->typeLUT[0]->counts-1);
    free(samBlock->aux->models->typeLUT[0]);
    free(samBlock->aux->models->typeLUT);
    free(samBlock->aux->models->typeRAW[0]->counts-1);
    free(samBlock->aux->models->typeRAW[0]);
    free(samBlock->aux->models->typeRAW);
    free(samBlock->aux->models->descBytes[0]->counts-1);
    free(samBlock->aux->models->descBytes[0]);
    free(samBlock->aux->models->descBytes[1]->counts-1);
    free(samBlock->aux->models->descBytes[1]);
    free(samBlock->aux->models->descBytes);
    free(samBlock->aux->models->iidBytes[0]->counts-1);
    free(samBlock->aux->models->iidBytes[0]);
    free(samBlock->aux->models->iidBytes);
    free(samBlock->aux->models->most_common_values[0]->counts-1);
    free(samBlock->aux->models->most_common_values[0]);
    free(samBlock->aux->models->most_common_values);
    free(samBlock->aux->models->most_common_flag[0]->counts-1);
    free(samBlock->aux->models->most_common_flag[0]);
    free(samBlock->aux->models->most_common_flag);
    free(samBlock->aux->models->most_common_list[0]->counts-1);
    free(samBlock->aux->models->most_common_list[0]);
    free(samBlock->aux->models->most_common_list);
    free(samBlock->aux->models->sign_integers[0]->counts-1);
    free(samBlock->aux->models->sign_integers[0]);
    free(samBlock->aux->models->sign_integers);
    for (int i = 0; i < 4; i++) {
        free(samBlock->aux->models->integers[i]->counts-1);
        free(samBlock->aux->models->integers[i]);
    }
    free(samBlock->aux->models->integers);
    for (int i = 0; i < MAXLUT + 2; i++){
        free(samBlock->aux->models->aux_TagType[i]->counts-1);
        free(samBlock->aux->models->aux_TagType[i]);
        free(samBlock->aux->models->most_common_values_wContext[i]->counts-1);
        free(samBlock->aux->models->most_common_values_wContext[i]);
        free(samBlock->aux->models->sign_integers_wContext[i]->counts-1);
        free(samBlock->aux->models->sign_integers_wContext[i]);
        free(samBlock->aux->models->iidBytes_wContext[i]->counts-1);
        free(samBlock->aux->models->iidBytes_wContext[i]);
        for (int j = 0; j < 4; j++){
            free(samBlock->aux->models->integers_wContext[i*4 + j]->counts-1);
            free(samBlock->aux->models->integers_wContext[i*4 + j]);
        }
        for (int j = 0; j < 2; j++){
            free(samBlock->aux->models->descBytes_wContext[i*2 + j]->counts-1);
            free(samBlock->aux->models->descBytes_wContext[i*2 + j]);
        }
    }
    free(samBlock->aux->models->aux_TagType);
    free(samBlock->aux->models->most_common_values_wContext);
    free(samBlock->aux->models->sign_integers_wContext);
    free(samBlock->aux->models->iidBytes_wContext);
    free(samBlock->aux->models->integers_wContext);
    free(samBlock->aux->models->descBytes_wContext);
    free(samBlock->aux->models->firstAux_TagType[0]->counts-1);
    free(samBlock->aux->models->firstAux_TagType[0]);
    free(samBlock->aux->models->firstAux_TagType);
    free(samBlock->aux->models);
    free(samBlock->aux);
    for (int i = 0; i < MAX_LINES_PER_BLOCK; i++)
        free(samBlock->rnames->rnames[i]);
    free(samBlock->rnames->rnames);
    free(samBlock->rnames->models->same_ref[0]->counts-1);
    free(samBlock->rnames->models->same_ref[0]);
    free(samBlock->rnames->models->same_ref);
    for (int i = 0; i < 256; i++){
        free(samBlock->rnames->models->rname[i]->counts-1);
        free(samBlock->rnames->models->rname[i]);
    }
    free(samBlock->rnames->models->rname);
    free(samBlock->rnames->models);
    free(samBlock->rnames);
    for (int i = 0; i < 256; i++){
        free(samBlock->mapq->models->mapq[i]->counts-1);
        free(samBlock->mapq->models->mapq[i]);
    }
    free(samBlock->mapq->models->mapq);
    free(samBlock->mapq->models);
    free(samBlock->mapq->mapq);
    free(samBlock->mapq);
    free(samBlock->rnext->rnext[0]);
    free(samBlock->rnext->rnext);
    free(samBlock->rnext->models->same_ref[0]->counts-1);
    free(samBlock->rnext->models->same_ref[0]);
    free(samBlock->rnext->models->same_ref);
    for (int i = 0; i < 256; i++){
        free(samBlock->rnext->models->rnext[i]->counts-1);
        free(samBlock->rnext->models->rnext[i]);
    }
    free(samBlock->rnext->models->rnext);
    free(samBlock->rnext->models);
    free(samBlock->rnext);
    free(samBlock->pnext->pnext);
    free(samBlock->pnext->models->zero[0]->counts-1);
    free(samBlock->pnext->models->zero[0]);
    free(samBlock->pnext->models->zero);
    free(samBlock->pnext->models->sign[0]->counts-1);
    free(samBlock->pnext->models->sign[0]);
    free(samBlock->pnext->models->sign);
    free(samBlock->pnext->models->assumption[0]->counts-1);
    free(samBlock->pnext->models->assumption[0]);
    free(samBlock->pnext->models->assumption);
    for (int i = 0; i < 4; i++){
        free(samBlock->pnext->models->raw_pnext[i]->counts-1);
        free(samBlock->pnext->models->raw_pnext[i]);
        free(samBlock->pnext->models->diff_pnext[i]->counts-1);
        free(samBlock->pnext->models->diff_pnext[i]);
    }
    free(samBlock->pnext->models->raw_pnext);
    free(samBlock->pnext->models->diff_pnext);
    free(samBlock->pnext->models);
    free(samBlock->pnext);
    free(samBlock->tlen->tlen);
    free(samBlock->tlen->models->sign[0]->counts-1);
    free(samBlock->tlen->models->sign[0]);
    free(samBlock->tlen->models->sign);
    free(samBlock->tlen->models->zero[0]->counts-1);
    free(samBlock->tlen->models->zero[0]);
    free(samBlock->tlen->models->zero);
    for (int i = 0; i < 4; i++){
        free(samBlock->tlen->models->tlen[i]->counts-1);
        free(samBlock->tlen->models->tlen[i]);
    }
    free(samBlock->tlen->models->tlen);
    free(samBlock->tlen->models);
    free(samBlock->tlen);
    free(samBlock->QVs);
    free(samBlock);
 
}

void* compress(void *thread_info){
    char prefix[2]; prefix[0] = 'F'; prefix[1] = 'F';
    uint64_t context[25][6];

    for(int i = 0; i < 25; i++) {
        for(int j = 0; j < 6; j++){
            if (j == 5) { context[i][j] = 5; }
            else { context[i][j] = 1; }
        }
    }

    uint64_t compress_file_size = 0;
    uint64_t compress_file_size1 = 0;
    clock_t begin;
    clock_t ticks;
    
    unsigned long long lineCtr = 0;
    
    printf("Compressing...\n");
    begin = clock();
    
    struct compressor_info_t info = *((struct compressor_info_t *)thread_info);
    
    struct qv_options_t opts = *(info.qv_opts);
    char buffer[1024];
    // Allocs the Arithmetic and the I/O stream
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);
    Arithmetic_stream as1 = alloc_arithmetic_stream(info.mode, info.frefcom);    

    // Allocs the different blocks and all the models for the Arithmetic
    sam_block samBlock = alloc_sam_models(as, info.fsam, info.fref, &opts, info.mode);
    
    //idoia
    create_most_common_list(samBlock);
    compress_most_common_list(as, samBlock->aux);
    //idoia
    
    if (info.lossiness == LOSSY) {
        compress_int(as, samBlock->codebook_model, LOSSY);
        initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    }
    else
        compress_int(as, samBlock->codebook_model, LOSSLESS);
   
    int polyploidy_ = 2;
    int qualityValueMax_ = 41;
    int qualityValueMin_ = 0;
    int qualityValueOffset_ = 33;


    calq::CQFile cqFile("quality_values_calq", calq::File::MODE_WRITE);
    std::cout << samBlock->block_length << std::endl;
    cqFile.writeHeader(10000);

    printf("start line compression\n"); 
    std::vector<calq::SAMRecord> samRecords;
    while (compress_line(as, samBlock, info.funmapped, info.lossiness, info.calqmode, as1, context, prefix, info.fsinchr, samRecords)) {
        ++lineCtr;
        if (lineCtr % 10000 == 0) {
          std::cout << "Starting qualEncoder with " << lineCtr << " lines." << std::endl;
          lineCtr = 0;
          if (info.calqmode){
            calq::QualEncoder qualEncoder(polyploidy_, qualityValueMax_, qualityValueMin_, qualityValueOffset_);
            std::cout << "initialized qual encoder" << std::endl;
            int i = 0;
            for (auto const &samRecord : samRecords) {
                    std::cout << i << "/" << samRecords.size() << " ";
                    qualEncoder.addMappedRecordToBlock(samRecord);
                    i++;
            }
            std::cout << "Sanity check" << std::endl;
            qualEncoder.finishBlock();
            qualEncoder.writeBlock(&cqFile);
          }
          printf("[cbc] compressed %llu lines\n", lineCtr);
          samRecords.clear();
        }
    }
    if (info.calqmode){
      calq::QualEncoder qualEncoder(polyploidy_, qualityValueMax_, qualityValueMin_, qualityValueOffset_);
      for (auto const &samRecord : samRecords) {
          qualEncoder.addMappedRecordToBlock(samRecord);
      }
      qualEncoder.finishBlock();
      qualEncoder.writeBlock(&cqFile);
    }
    printf("[cbc] compressed %llu lines\n", lineCtr);

    samRecords.clear();
    printf("done compressing lines\n"); 
    // Check if we are in the last block
    compress_rname(as, samBlock->rnames->models, "\n");
  
    //end the compression
    compress_file_size = encoder_last_step(as);
    compress_file_size1 = encoder_last_step(as1);    
 
    printf("Final Size: %lu\n", compress_file_size);
    printf("Final Size (Compressed Reference File): %lu\n", compress_file_size1);

    ticks = clock() - begin;
    
    printf("Compression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);

    free_os_stream(as1->ios);
    free(as1);
    free_os_stream(as->ios);
    free(as);
    free_sam_block_compress(samBlock);
    //pthread_exit(NULL);

    return NULL;
}


void* decompress(void *thread_info){
    uint64_t n = 0;
    clock_t begin = clock();
    clock_t ticks;
    

    printf("Decompression started \n");
    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    printf("Re-Initialized compressor info\n");
    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
    calq::CQFile fcq_(info->cq_name, calq::CQFile::MODE_READ);
    size_t blockSize = 0;  
    calq::File qualFile_("Calq_Output.txt", calq::File::MODE_WRITE);

    if (info->calqmode){
        printf("Reading cq file Info of %s\n", info->cq_name);
        fcq_.readHeader(&blockSize);
    }

    sam_block samBlock = alloc_sam_models(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    decompress_most_common_list(as, samBlock->aux);
    
    info->lossiness = decompress_int(as, samBlock->codebook_model);
    
    // Start the decompression
    // initialize the QV model
    if (info->lossiness == LOSSY) {
        initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    }
    // Decompress the blocks
    std::vector<calq::SAMRecord> samRecords;
    while(decompress_line(as, samBlock, info->lossiness, info->calqmode, samRecords)){
        n++;

        if (n % blockSize == 0 && info->calqmode) {
            calq::QualDecoder qualDecoder;
            qualDecoder.readBlock(&fcq_);
            for (auto const &samRecord : samRecords) {
                qualDecoder.decodeMappedRecordFromBlock(samRecord, &qualFile_);
            }
        }
    }
    if (info->calqmode) {
        calq::QualDecoder qualDecoder;
       qualDecoder.readBlock(&fcq_);

        printf("decoding..:");
        for (auto const &samRecord : samRecords) {
            qualDecoder.decodeMappedRecordFromBlock(samRecord, &qualFile_);
        }
    }
    
    n += samBlock->block_length;
    free_sam_block_compress(samBlock);
    ticks = clock() - begin;
    printf("Decompression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    return NULL;
}
