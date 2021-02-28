
import sys
import logging
from ksrates.utils import init_logging
from ksrates.cluster_anchor_ks import cluster_anchor_ks
from ksrates.exp_log_mixture import exp_log_mixture
from ksrates.lognormal_mixture import lognormal_mixture
import ksrates.fc_configfile as fcConf

def paralogs_analyses_methods(config_file, paranome_table, anchors_table, correction_table, anchorpoints, multiplicons, segments, list_elements, multiplicon_pairs):
  # INPUT
  config = fcConf.Configuration(config_file)
  logging.basicConfig(format='%(levelname)s\t%(message)s', level=config.get_logging_level(), stream=sys.stdout)

  paranome = config.get_paranome()
  colinearity = config.get_colinearity()
  extra_paralogs_analyses_methods = config.get_extra_paralogs_analyses_methods()
  
  if paranome and not colinearity:
    # Only exp-log mixture model by default
    exp_log_mixture(config_file, paranome_table, correction_table)
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      # Lognormal mixture model on paranome
      lognormal_mixture(config_file, paranome_table, anchors_table, correction_table) 
      
  if colinearity and not paranome:
    # Only anchor clustering by default
    cluster_anchor_ks(config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      # Lognormal mixture model on anchors
      lognormal_mixture(config_file, paranome_table, anchors_table, correction_table)

  if colinearity and paranome:
    # Only anchor clustering by default
    cluster_anchor_ks(config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      # Exp-log mixture model on paranome
      exp_log_mixture(config_file, paranome_table, correction_table)
      logging.info(f"\n")
      # Lognormal mixture model on both
      lognormal_mixture(config_file, paranome_table, anchors_table, correction_table)
