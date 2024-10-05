"use client"
import React from 'react';
import { createViewState, JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view';

export const GenomeAnalyser = () => {
    const assembly = {
        name: 'GRCh38',
        sequence: {
            type: 'ReferenceSequenceTrack',
            trackId: 'GRCh38-ReferenceSequenceTrack',
            adapter: {
                type: 'BgzipFastaAdapter',
                fastaLocation: {
                    uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz',
                },
                faiLocation: {
                    uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz.fai',
                },
                gziLocation: {
                    uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz.gzi',
                },
            },
        },
        aliases: ['hg38'],
        refNameAliases: {
            adapter: {
                type: 'RefNameAliasAdapter',
                location: {
                    uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt',
                },
            },
        },
    };


    const tracks = [
        {
            type: 'VariantTrack',
            trackId: '1000genomes_variants',
            name: '1000 Genomes Variants',
            assemblyNames: ['GRCh38'],
            adapter: {
                type: 'VcfTabixAdapter',
                vcfGzLocation: {
                    uri: 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz',
                },
                index: {
                    location: {
                        uri: 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi',
                    },
                },
            },
        },
    ];


    const state = createViewState({
        assembly,
        tracks,
        location: '17:43,044,295..43,125,482', // BRCA1 location
    });

    // Highlight BRCA1 and BRCA2 regions
    const highlightedRegions = [
        {
            assemblyName: 'GRCh38',
            refName: '17',
            start: 43044295,
            end: 43125482,
            reversed: false,
        },
        {
            assemblyName: 'GRCh38',
            refName: '13',
            start: 32315474,
            end: 32400266,
            reversed: false,
        },
    ];
    state.session.view.setHighlight(highlightedRegions);

    return (
        <div>
            <JBrowseLinearGenomeView viewState={state} />
        </div>
    );
};