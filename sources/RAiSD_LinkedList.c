/*  
 *  RAiSD, Raised Accuracy in Sweep Detection
 *
 *  Copyright January 2017 by Nikolaos Alachiotis and Pavlos Pavlidis
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an email to
 *  Nikolaos Alachiotis (n.alachiotis@gmail.com)
 *  Pavlos Pavlidis (pavlidisp@gmail.com)  
 *  
 */

#include "RAiSD.h"

RSDLinkedListNode_t * RSDLinkedList_new (double pos, char * snp)
{
	RSDLinkedListNode_t * newNode = NULL;

	newNode = (RSDLinkedListNode_t *)malloc(sizeof(RSDLinkedListNode_t));
	assert(newNode!=NULL);

	newNode->prv = NULL;
	newNode->nxt = NULL;
	newNode->snpPosition = pos;
	
	int snpSize = strlen(snp)+1;
	newNode->snpData = (char *)malloc(sizeof(char)*snpSize);
	assert(newNode->snpData!=NULL);

	strncpy(newNode->snpData, snp, strlen(snp));
	newNode->snpData[strlen(snp)] = '\0';

	return newNode;
}

// redo simply linked
RSDLinkedListNode_t * RSDLinkedList_addNode (RSDLinkedListNode_t * listHead, double pos, char * snp)
{
	assert(snp!=NULL);
	assert(pos>=0.0);

	RSDLinkedListNode_t * newNode = RSDLinkedList_new (pos, snp);
	assert(newNode!=NULL);
	
	if(listHead==NULL) // empty list
	{
		return newNode;	
	}
	else
	{
		RSDLinkedListNode_t * curNode = listHead;
		assert(curNode!=NULL);

		RSDLinkedListNode_t * lstNode = NULL;

		double curPos = curNode->snpPosition;

		// Handle case of new head
		if(pos<curPos)
		{
			newNode->nxt = curNode;
			curNode->prv = newNode;

			return newNode; // this becomes the new list head
		}
		else
		{
			int newNodeAdded = 0;

			while(curNode!=NULL)
			{
				curPos = curNode->snpPosition;

				if(curNode->nxt==NULL)
					lstNode = curNode;
			
				if(curPos>pos)
				{
					newNode->nxt = curNode;
					newNode->prv = curNode->prv;
					newNode->prv->nxt = newNode;
					curNode->prv = newNode;
					
					newNodeAdded = 1;

					break;
				}

				curNode = curNode->nxt;
			}

			if(newNodeAdded==0)
			{
				assert(lstNode!=NULL);

				lstNode->nxt = newNode;
				newNode->prv = lstNode;

				lstNode = NULL;
			}

			return listHead;
		}		
	}
	
	return listHead;
}

int RSDLinkedList_getSize (RSDLinkedListNode_t * listHead)
{
	assert(listHead!=NULL);
	int sz = 0;
	
	RSDLinkedListNode_t * curNode = listHead;

	while(curNode!=NULL)
	{
		sz++;
		curNode = curNode->nxt;
	}

	return sz;
}

void RSDLinkedList_appendToFile (RSDLinkedListNode_t * listHead, FILE * fp)
{
	assert(listHead!=NULL);
	assert(fp!=NULL);

	RSDLinkedListNode_t * curNode = listHead;

	while(curNode!=NULL)
	{
		fprintf(fp, "%s\n", curNode->snpData);
		curNode = curNode->nxt;
	}
	fflush(fp);	
}

RSDLinkedListNode_t * RSDLinkedList_free(RSDLinkedListNode_t * listHead)
{
	assert(listHead!=NULL);

	RSDLinkedListNode_t * curNode = listHead;
	RSDLinkedListNode_t * nxtNode = curNode->nxt;

	while(curNode!=NULL)
	{
		if(curNode->snpData!=NULL)
		{
			free(curNode->snpData);
			curNode->snpData=NULL;
		}
		
		nxtNode = curNode->nxt;

		free(curNode);
		
		curNode = nxtNode;
	}

	return NULL;
}














