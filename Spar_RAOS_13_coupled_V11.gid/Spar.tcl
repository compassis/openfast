
proc TdynTcl_InitiateProblem { } {


    configure_analysis Solve_Dif_Rad 1
    #  TdynTcl_Message "Solve_dif_rad set to 0!" notice
	###### Mass matrix of spar buoy OC3 ######
	# TdynTcl_Add_Mass_Matrix 1 [list 0.0 0.0 0.0] [list \
	# 7.4663E+06		0				0				0				-6.7134E+08		0 \
	# 0	 				7.4663E+06		0				6.7134E+08		0				0 \
	# 0			  		0				7.4663E+06		0				0				0 \
	# 0	 				6.7134E+08		0				6.4593E+10		0				0 \
    # -6.7134E+08		0				0				0				6.4593E+10		0 \
	# 0	 				0				0				0				0				1.6420E+08 ]
	
	
	###### Mass matrix of spar buoy OC3 ######
	# TdynTcl_Add_Mass_Matrix 1 [list 0.0 0.0 -89.9155] [list \
	# 7.46633E+06			0				0				0				0				0 \
	# 0	 				7.46633E+06		0				0				0				0 \
	# 0			  		0				7.46633E+06		0				0				0 \
	# 0	 				0				0				4.22923E-09		0				0 \
    # 0					0				0				0				4.22923E-09		0 \
	# 0	 				0				0				0				0				1.6423E-08 ]
	   
	###### Additional Damping of spar buoy OC3 ######
	TdynTcl_Add_Damping_Matrix 1 [list 0.0 0.0 0.0] [list \
	 100000.0			0		 	0		  0	  		 0			 0 \
			0	 100000.0	       	0		  0	  		 0			 0 \
			0			0	 130000.0		  0	  	     0			 0 \
			0	 		0			0		  0	  	     0			 0 \
			0			0			0	  	  0			 0		     0 \
			0	 		0			0		  0	  		 0	   1.3E+07 ]
			 
	
	###### Mass matrix of wind turbine ######
	# TdynTcl_Add_Mass_Matrix 1 [list 0.0 0.0 0.0] [list \
	# 6.9746E+05			0			 0			 0		4.430E+07				 0 \
			 # 0	6.975E+05			 0	 -4.43E+07				0		6.6000E+06 \
			 # 0			0	6.9746E+05			 0	  -6.6000E+06			     0 \
			 # 0	-4.43E+07			 0	3499.0E+06				0	   			 0 \
	  # 4.43E+07			0	 -6.60E+06			 0	   3560.0E+06			     0 \
			 # 0	6.600E+06			 0	 -5.133E+08				0		1.0117E+08 ]

	
	###### Damping matrix of wind turbine ######
	# TdynTcl_Add_Damping_Matrix 1 [list 0.0 0.0 0.0] [list \
    # 0.04E+06            0        -0.01E+04        -0.25E+06        4.00E+06        0.08E+06 \
		   # 0     		 0                0        -0.11E+06       -0.18E+06       -0.05E+06 \
	# -0.01E+06            0         		  0        -0.04E+06       -0.92E+06       -0.33E+06 \
	 # 0.27E+06    -0.10E+06                0        16.17E+06       50.30E+06       13.88E+06 \
	 # 3.42E+06     0.06E+06        -1.00E+06       -23.92E+06       400.1E+06       59.01E+06 \
	 # 0.05E+06    -0.02E+06         0.22E+06       11.08E+06      -52.60E+06       102.2E+06 ]
   
	###### Stiffness matrix of wind turbine ######
    # TdynTcl_Add_Stiffness_Matrix 1 [list 0.0 0.0 0.0] [list \
		   # 0        0        0        0.3E+06        0.2E+06                0 \
		   # 0        0        0       -0.1E+06        0.3E+06        -0.07E+06 \
		   # 0        0        0       -0.3E+06       -0.4E+06                0 \
		   # 0        0        0        8.5E+06      -22.4E+06         59.7E+06 \
		   # 0        0        0       26.8E+06       28.9E+06         -4.1E+06 \
		   # 0        0        0       -1.2E+06        1.1E+06         -4.8E+06 ]
	

	###### Stiffness matrix for mooring linear analysis ######
	TdynTcl_Add_Stiffness_Matrix 1 [list 0.0 0.0 0.0] [list \
	 4.112E+04              0         0              0     -2.821E+06            0 \
	         0      4.112E+04         0      2.821E+06              0            0 \
	         0              0 1.194E+04              0              0            0 \
	         0      2.816E+06         0      3.111E+08              0            0 \
	-2.816E+06              0         0              0      3.111E+08            0 \
	         0              0         0              0              0    11.56E+06 ]
	
	TdynTcl_Message "TdynTcl_DefineBodyData finished!!!" notice
}