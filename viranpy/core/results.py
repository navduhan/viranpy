#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Result classes for ViRAnPy pipeline.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional
from pathlib import Path
import json


@dataclass
class AnnotationResult:
    """
    Container for annotation results.
    
    This class holds the results from various annotation steps
    such as gene prediction, tRNA detection, etc.
    """
    
    annotator_name: str
    input_file: str
    output_files: List[str] = field(default_factory=list)
    annotations: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    success: bool = True
    error_message: Optional[str] = None
    
    def add_annotation(self, key: str, value: Any) -> None:
        """
        Add an annotation to the results.
        
        Args:
            key: Annotation key
            value: Annotation value
        """
        self.annotations[key] = value
    
    def add_output_file(self, file_path: str) -> None:
        """
        Add an output file to the results.
        
        Args:
            file_path: Path to output file
        """
        self.output_files.append(file_path)
    
    def add_metadata(self, key: str, value: Any) -> None:
        """
        Add metadata to the results.
        
        Args:
            key: Metadata key
            value: Metadata value
        """
        self.metadata[key] = value
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert results to dictionary.
        
        Returns:
            Dictionary representation of results
        """
        return {
            'annotator_name': self.annotator_name,
            'input_file': self.input_file,
            'output_files': self.output_files,
            'annotations': self.annotations,
            'metadata': self.metadata,
            'success': self.success,
            'error_message': self.error_message,
        }
    
    def save(self, output_file: str) -> None:
        """
        Save results to a JSON file.
        
        Args:
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, input_file: str) -> 'AnnotationResult':
        """
        Load results from a JSON file.
        
        Args:
            input_file: Input file path
            
        Returns:
            AnnotationResult instance
        """
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        return cls(**data)


@dataclass
class PredictionResult:
    """
    Container for prediction results.
    
    This class holds the results from various prediction steps
    such as protein function prediction, etc.
    """
    
    predictor_name: str
    input_file: str
    output_files: List[str] = field(default_factory=list)
    predictions: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    success: bool = True
    error_message: Optional[str] = None
    
    def add_prediction(self, key: str, value: Any) -> None:
        """
        Add a prediction to the results.
        
        Args:
            key: Prediction key
            value: Prediction value
        """
        self.predictions[key] = value
    
    def add_output_file(self, file_path: str) -> None:
        """
        Add an output file to the results.
        
        Args:
            file_path: Path to output file
        """
        self.output_files.append(file_path)
    
    def add_metadata(self, key: str, value: Any) -> None:
        """
        Add metadata to the results.
        
        Args:
            key: Metadata key
            value: Metadata value
        """
        self.metadata[key] = value
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert results to dictionary.
        
        Returns:
            Dictionary representation of results
        """
        return {
            'predictor_name': self.predictor_name,
            'input_file': self.input_file,
            'output_files': self.output_files,
            'predictions': self.predictions,
            'metadata': self.metadata,
            'success': self.success,
            'error_message': self.error_message,
        }
    
    def save(self, output_file: str) -> None:
        """
        Save results to a JSON file.
        
        Args:
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, input_file: str) -> 'PredictionResult':
        """
        Load results from a JSON file.
        
        Args:
            input_file: Input file path
            
        Returns:
            PredictionResult instance
        """
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        return cls(**data)


@dataclass
class PipelineResult:
    """
    Container for complete pipeline results.
    
    This class holds the results from all pipeline steps
    and provides methods for summary and export.
    """
    
    config: Dict[str, Any]
    annotation_results: List[AnnotationResult] = field(default_factory=list)
    prediction_results: List[PredictionResult] = field(default_factory=list)
    output_files: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    success: bool = True
    error_message: Optional[str] = None
    
    def add_annotation_result(self, result: AnnotationResult) -> None:
        """
        Add an annotation result.
        
        Args:
            result: Annotation result to add
        """
        self.annotation_results.append(result)
    
    def add_prediction_result(self, result: PredictionResult) -> None:
        """
        Add a prediction result.
        
        Args:
            result: Prediction result to add
        """
        self.prediction_results.append(result)
    
    def add_output_file(self, file_path: str) -> None:
        """
        Add an output file.
        
        Args:
            file_path: Path to output file
        """
        self.output_files.append(file_path)
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of pipeline results.
        
        Returns:
            Dictionary containing summary information
        """
        return {
            'success': self.success,
            'total_annotations': len(self.annotation_results),
            'total_predictions': len(self.prediction_results),
            'total_output_files': len(self.output_files),
            'successful_annotations': sum(1 for r in self.annotation_results if r.success),
            'successful_predictions': sum(1 for r in self.prediction_results if r.success),
            'error_message': self.error_message,
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert results to dictionary.
        
        Returns:
            Dictionary representation of results
        """
        return {
            'config': self.config,
            'annotation_results': [r.to_dict() for r in self.annotation_results],
            'prediction_results': [r.to_dict() for r in self.prediction_results],
            'output_files': self.output_files,
            'metadata': self.metadata,
            'success': self.success,
            'error_message': self.error_message,
            'summary': self.get_summary(),
        }
    
    def save(self, output_file: str) -> None:
        """
        Save results to a JSON file.
        
        Args:
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, input_file: str) -> 'PipelineResult':
        """
        Load results from a JSON file.
        
        Args:
            input_file: Input file path
            
        Returns:
            PipelineResult instance
        """
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        # Reconstruct annotation results
        annotation_results = [
            AnnotationResult(**r) for r in data.get('annotation_results', [])
        ]
        
        # Reconstruct prediction results
        prediction_results = [
            PredictionResult(**r) for r in data.get('prediction_results', [])
        ]
        
        return cls(
            config=data['config'],
            annotation_results=annotation_results,
            prediction_results=prediction_results,
            output_files=data.get('output_files', []),
            metadata=data.get('metadata', {}),
            success=data.get('success', True),
            error_message=data.get('error_message'),
        ) 