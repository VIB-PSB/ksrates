
import sys
import logging
from ksrates.cluster_anchor_ks import cluster_anchor_ks
from ksrates.exp_log_mixture import exp_log_mixture
from ksrates.lognormal_mixture import lognormal_mixture
import ksrates.fc_configfile as fcConf

def paralogs_analyses_methods(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table,
                  correction_table, anchorpoints, multiplicons, segments, list_elements, multiplicon_pairs):
  # INPUT
  config = fcConf.Configuration(config_file, expert_config_file)
  logging.basicConfig(format='%(levelname)s\t%(message)s', level=config.get_logging_level(), stream=sys.stdout)

  paranome = config.get_paranome()
  colinearity = config.get_colinearity()
  reciprocal_retention_analysis = config.get_reciprocal_retention()
  extra_paralogs_analyses_methods = config.get_extra_paralogs_analyses_methods()
  
  # 1) ONLY WHOLE-PARANOME
  if paranome and not colinearity and not reciprocal_retention_analysis:
    # Default: perform only exponential-lognormal mixture modeling (ELMM) on paranome
    exp_log_mixture(config_file, expert_config_file, paranome_table, correction_table)

    # If additional methods are asked: perform lognormal mixture modeling (LMM) on paranome
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table) 
    

  # 2) ONLY COLLINEARITY 
  elif not paranome and colinearity and not reciprocal_retention_analysis:
    # Default: perform only anchor clustering
    cluster_anchor_ks(config_file, expert_config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)

    # If additional methods are asked: perform LMM on anchors
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)


  # 3) ONLY REC.RET.
  elif not paranome and not colinearity and reciprocal_retention_analysis:
    # Default: perform LMM on rec_ret paralogs
    lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)
    # There are no additional methods for rec_ret paralogs

  
  # 4) COLLINEARITY & REC.RET.  
  elif not paranome and colinearity and reciprocal_retention_analysis:
    # Default: perform 1) anchor clustering and 2) LMM on rec_ret paralogs (see "if" block below)
    cluster_anchor_ks(config_file, expert_config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)
    logging.info(f"\n")
    if not extra_paralogs_analyses_methods: # Note: overwrite anchors to be ignored
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table, ignore_anchors=True)

    # If additional methods are asked: perform LMM on rec_ret paralogs (like in default) plus anchors
    elif extra_paralogs_analyses_methods:
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)


  # 5) COLLINEARITY & PARANOME
  elif paranome and colinearity and not reciprocal_retention_analysis:
    # Default: perform only anchor clustering
    cluster_anchor_ks(config_file, expert_config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)

    # If additional methods are asked: perform ELMM on paranome, and LMM on paranome and anchors
    if extra_paralogs_analyses_methods:
      logging.info(f"\n")
      # Exp-log mixture model on paranome
      exp_log_mixture(config_file, expert_config_file, paranome_table, correction_table)
      logging.info(f"\n")
      # Lognormal mixture model on paranome and anchors
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)


  # 6) REC.RET. & PARANOME
  elif paranome and not colinearity and reciprocal_retention_analysis:
      # Default: perform ELMM on paranome and LMM on rec_ret paralogs
      exp_log_mixture(config_file, expert_config_file, paranome_table, correction_table)
      logging.info(f"\n")
      if not extra_paralogs_analyses_methods: # Note: overwrite paranome to be ignored
        lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table, ignore_paranome=True)

      # If additional methods are asked: perform LMM on rec_ret paralogs (like in default) plus paranome
      elif extra_paralogs_analyses_methods:
        lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)


  # 7) COLLINEARITY & REC.RET. & PARANOME
  elif paranome and colinearity and reciprocal_retention_analysis:
    # Default: perform anchor clustering and LMM on rec_ret paralogs
    cluster_anchor_ks(config_file, expert_config_file, correction_table, anchorpoints, multiplicons, segments, list_elements, anchors_table, multiplicon_pairs)
    logging.info(f"\n")
    if not extra_paralogs_analyses_methods: # Note: overwrite paranome and anchors to be ignored
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table, ignore_anchors=True, ignore_paranome=True)

    # If additional methods are asked: perform ELMM on paranome, and LMM on rec_ret paralogs (like in default), plus paranome and anchors
    elif extra_paralogs_analyses_methods:
      # Exp-log mixture model on paranome
      exp_log_mixture(config_file, expert_config_file, paranome_table, correction_table)
      logging.info(f"\n")
      # Lognormal mixture model on anchors, paranome and rec_ret
      lognormal_mixture(config_file, expert_config_file, paranome_table, anchors_table, reciprocal_retention_table, correction_table)
