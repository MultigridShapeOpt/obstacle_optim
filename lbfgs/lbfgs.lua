local deque = require "deque"

local init_grid_functions= function(self, AppSpace)
	if AppSpace ~= nil then 
	print("Given Approximation Space, initialize GridFunction vectors")
		self:createGFVector(AppSpace)
		print("Create individual gfs")
		self.q_tmp=GridFunction(AppSpace);
		self.q_tmp:set(0.0);
		self.x_old=GridFunction(AppSpace);
		self.x_old:set(0.0);
		self.g_old=GridFunction(AppSpace);
		self.g_old:set(0.0);
		self.inner_gfs[1]=GridFunction(AppSpace);self.inner_gfs[1]:set(0.0);
		self.inner_gfs[2]=GridFunction(AppSpace);self.inner_gfs[2]:set(0.0);
		self.inner_gfs[3]=GridFunction(AppSpace);self.inner_gfs[3]:set(0.0);
	else 
		print("BFGS: Approximation Space is null, GridFunctions not initialized")
	end
	
end

local createGFVector = function(self, AppSpace)
	--self.q_tmp=GridFunction(AppSpace);self.x_old=GridFunction(AppSpace);self.g_old=GridFunction(AppSpace);
	--self.q_tmp:set(0.0);self.x_old:set(0.0);self.g_old:set(0.0)
	for i=1, self.memory do
		self.y[i]=GridFunction(AppSpace);self.y[i]:set(0.0);
		self.s[i]=GridFunction(AppSpace);self.s[i]:set(0.0);
	end
end

local reset = function(self)
	self.iter=0
	self.indeces:clear()
end

local print_status = function(self)
	print("Size of y "..#self.y)
	print("Size of s "..#self.s)
	print("Value of iter "..self.iter)
	print("Value of memory "..self.memory)
	print("Size of rho "..#self.rho)
	print("Size of alpha "..#self.valpha)
end


local update_internal = function(self,x_k,g_k)
	if self.iter>0 and self.memory>0 then
	--inner_gfs[1] --> q_tmp
	--inner_gfs[2] --> x_old
	--inner_gfs[3] --> g_old
		--print("DO SOMETHING CASE")
		--print("WE PUSH WITH ITER "..self.iter)
		VecScaleAssign(self.inner_gfs[1],1.0,g_k)--self.inner_gfs[1]=g_k:clone()
		--update and store x
		VecScaleAdd2(self.inner_gfs[2],1.0,x_k,-1.0,self.inner_gfs[2])
		i=self.iter%self.memory
		if i==0 then 
			i=self.memory 
		end
		--print("The value of iter is "..self.iter)
		--print("The value of i is: "..i)
		VecScaleAssign(self.s[i],1.0,self.inner_gfs[2])--self.s[i]=(self.inner_gfs[2]):clone()
		self.indeces:push_right(i)
		--size control
		if self.indeces:length() > self.memory then 
			self.indeces:pop_left()
		end
		--store old value
		VecScaleAssign(self.inner_gfs[2],1.0,x_k)--self.inner_gfs[2]=x_k:clone()
		
		--update and store g
		VecScaleAdd2(self.inner_gfs[3],1.0,g_k,-1.0,self.inner_gfs[3])
		VecScaleAssign(self.y[i],1.0,self.inner_gfs[3])--self.y[i]=(self.inner_gfs[3]):clone()
		VecScaleAssign(self.inner_gfs[3],1.0,g_k)--self.inner_gfs[3]=g_k:clone()

		--for rho, we need the most recent stores gridfuncs	
		--print("Dot product result is: "..self:call_dot_product(self.y[i],self.s[i]))
		self.rho[i]=1/self:call_dot_product(self.y[i],self.s[i])
		--print("rho is: "..self.rho[i])
		
		
		--cycle upwards
		for j in self.indeces:iter_right() do
			--print("DOWNWARDS CYCLE-STARTING ITERATION ")
			--print("The value of j is: "..j)
			self.valpha[j]=self.rho[j]*self:call_dot_product(self.inner_gfs[1],self.s[j])
			VecScaleAdd2(self.inner_gfs[1],1.0,self.inner_gfs[1],-self.valpha[j],self.y[j])
			--print("alpha is "..self.valpha[j])
		end
		--calculate gamma 
		--print("gamma numerator is: "..self:call_dot_product(self.y[i],self.s[i]))
		--print("gamma denominator is: "..self:call_dot_product(self.y[i],self.y[i]))
		local gamma=self:call_dot_product(self.y[i],self.s[i])/
					self:call_dot_product(self.y[i],self.y[i])
		--print("gamma is "..gamma)
		
		--calculate z (descent direction)
		VecScaleAssign(self.inner_gfs[1],gamma,self.inner_gfs[1])
		
		for j in self.indeces:iter_left() do
			--print("UPWARDS CYCLE-STARTING ITERATION "..idx)
			local beta = self:call_dot_product(self.inner_gfs[1],self.y[j])
			VecScaleAdd2(self.inner_gfs[1],1.0,self.inner_gfs[1],self.valpha[j]-beta,self.s[j])
			--print("beta is "..beta)
			--print("alpha[j]-beta is "..self.valpha[j]-beta)
		end	
		VecScaleAssign(g_k,-1.0,self.inner_gfs[1])
	elseif self.iter == 0 then
		print("INIT CASE")
		self.rho= {}--empty rho?
		VecScaleAssign(self.inner_gfs[2],1.0,x_k)--self.inner_gfs[2]=x_k:clone()
		VecScaleAssign(self.inner_gfs[3],1.0,g_k)--self.inner_gfs[3]=g_k:clone()
	else
		--print("DO NOTHING CASE")
	end
	self.iter=self.iter+1
end

local update = function(self,tmp,xown,gown,x_k,g_k)
	if self.iter>0 and self.memory>0 then
		print("DO SOMETHING CASE")
		print("WE PUSH WITH ITER "..self.iter)
		VecScaleAssign(tmp,1.0,g_k)--tmp=g_k:clone()
		--update and store x
		VecScaleAdd2(xown,1.0,x_k,-1.0,xown)
		i=self.iter%self.memory
		if i==0 then 
			i=self.memory 
		end
		print("The value of i is: "..i)
		VecScaleAssign(self.s[i],1.0,xown)--self.s[i]=xown:clone()
		--self.s[i]:assign(xown)
		self.indeces:push_right(i)
		--size control
		if self.indeces:length() > self.memory then 
			self.indeces:pop_left()
		end
		--store old value
		VecScaleAssign(xown,1.0,x_k)--xown=x_k:clone()
		--update and store g
		VecScaleAdd2(gown,1.0,g_k,-1.0,gown)
		VecScaleAssign(self.y[i],1.0,gown)--self.y[i]=gown:clone()
		VecScaleAssign(gown,1.0,g_k)--gown=g_k:clone()
		print("Dot product result for rho is: "..self:call_dot_product(self.y[i],self.s[i]))
		self.rho[i]=1/self:call_dot_product(self.y[i],self.s[i])
		print("rho is: "..self.rho[i])
		
		
		--cycle upwards
		for j in self.indeces:iter_right() do
			--print("DOWNWARDS CYCLE-STARTING ITERATION "..idx)
			print("The value of j is: "..j)
			self.valpha[j]=self.rho[j]*self:call_dot_product(tmp,self.s[j])
			VecScaleAdd2(tmp,1.0,tmp,-self.valpha[j],self.y[j])
			print("alpha is "..self.valpha[j])
		end
		--calculate gamma 					
		print("gamma numerator is: "..self:call_dot_product(self.y[i],self.s[i]))
		print("gamma denominator is: "..self:call_dot_product(self.y[i],self.y[i]))
		local gamma=self:call_dot_product(self.y[i],self.s[i])/
					self:call_dot_product(self.y[i],self.y[i])
					print("gamma is "..gamma)
		
		--calculate z (descent direction)
		VecScaleAssign(tmp,gamma,tmp)
		
		for j in self.indeces:iter_left() do
			local beta = self:call_dot_product(tmp,self.y[j])
			VecScaleAdd2(tmp,1.0,tmp,self.valpha[j]-beta,self.s[j])
			print("beta is "..beta)
			print("alpha[j]-beta is "..self.valpha[j]-beta)
		end	
		VecScaleAssign(g_k,-1.0,tmp)
	elseif self.iter == 0 then
		print("INIT CASE")
		self.rho= {}--empty rho?
		VecScaleAssign(xown,1.0,x_k)--xown= x_k:clone()
		VecScaleAssign(gown,1.0,g_k)--gown= g_k:clone()
		VecScaleAssign(tmp,1.0,g_k)--tmp= g_k:clone()--just for the disc
	else
		print("DO NOTHING CASE")
	end
	self.iter=self.iter+1
end

local update_all_external = function(self,tmp,xown,gown,x_k,g_k,vy,vs)
	if self.iter>0 and self.memory>0 then
		print("DO SOMETHING CASE")
		print("WE PUSH WITH ITER "..self.iter)
		VecScaleAssign(tmp,1.0,g_k)--tmp=g_k
		--update and store x
		VecScaleAdd2(xown,1.0,x_k,-1.0,xown)
		i=self.iter%self.memory
		if i==0 then 
			i=self.memory 
		end
		print("The value of i is: "..i)
		--vs[i]=xown:clone()
		VecScaleAssign(vs[i],1.0,xown)--vs[i]=xown
		--self.s[i]:assign(xown)
		self.indeces:push_right(i)
		--size control
		if self.indeces:length() > self.memory then 
			self.indeces:pop_left()
		end
		--store old value
		VecScaleAssign(xown,1.0,x_k)--xown=x_k
		--update and store g
		VecScaleAdd2(gown,1.0,g_k,-1.0,gown)
		--vy[i]=gown:clone()
		VecScaleAssign(vy[i],1.0,gown)--vy[i]=gown
		VecScaleAssign(gown,1.0,g_k)--gown=g_k
		print("Dot product result for rho is: "..self:call_dot_product(vy[i],vs[i]))
		self.rho[i]=1/self:call_dot_product(vy[i],vs[i])
		print("rho is: "..self.rho[i])
		
		
		--cycle upwards
		for j in self.indeces:iter_right() do
			--print("DOWNWARDS CYCLE-STARTING ITERATION "..idx)
			print("The value of j is: "..j)
			self.valpha[j]=self.rho[j]*self:call_dot_product(tmp,vs[j])
			VecScaleAdd2(tmp,1.0,tmp,-self.valpha[j],vy[j])
			print("alpha is "..self.valpha[j])
		end
		--calculate gamma 					
		print("gamma numerator is: "..self:call_dot_product(vy[i],vs[i]))
		print("gamma denominator is: "..self:call_dot_product(vy[i],vy[i]))
		local gamma=self:call_dot_product(vy[i],vs[i])/
					self:call_dot_product(vy[i],vy[i])
					print("gamma is "..gamma)
		
		--calculate z (descent direction)
		VecScaleAssign(tmp,gamma,tmp)
		
		for j in self.indeces:iter_left() do
			local beta = self:call_dot_product(tmp,vy[j])
			VecScaleAdd2(tmp,1.0,tmp,self.valpha[j]-beta,vs[j])
			print("beta is "..beta)
			print("alpha[j]-beta is "..self.valpha[j]-beta)
		end	
		VecScaleAssign(g_k,-1.0,tmp)
	elseif self.iter == 0 then
		print("INIT CASE")
		self.rho= {}--empty rho?
		VecScaleAssign(xown,1.0,x_k)--xown= x_k
		VecScaleAssign(gown,1.0,g_k)--gown= g_k
		VecScaleAssign(tmp,1.0,g_k)--tmp= g_k--just for the disc
	else
		print("DO NOTHING CASE")
	end
	self.iter=self.iter+1
end

local init = function(self,mmry,cmp,innerss,surfaces,scalar_product)
	self.memory=mmry
	self.valpha[mmry]=0
	self.rho[mmry]=0
	self.component=cmp 
	self.innerdomain=innerss
	self.boundarydomain=surfaces
	self.dot_product=scalar_product
end

local function_test = function(self, func)
	self.dot_product=func
	print("The result of the operation is:"..func(self.integration_order))
end

local call_dot_product= function(self,gf_a,gf_b)
	if self.boundarydomain ~= nil then 
		--print("Surface Integral")
		return self.dot_product(gf_a,gf_b,self.component,self.boundarydomain,self.innerdomain,self.integration_order)
	elseif self.boundarydomain == nil then
		--print("Volume Integral")
		return self.dot_product(gf_a,gf_b,self.component,self.innerdomain,self.integration_order)
	end
	print("Something else went wrong")
end

local print_vector_grid_functions = function(self)
	
	vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(self.component, "lbfgs_variable")
	for idx=1,self.memory do
		vtkWriter:print("bfgs_y", self.y[idx],idx,idx, false)
		vtkWriter:print("bfgs_s", self.s[idx],idx,idx, false)
	end
	--vtkWriter:print(ExtensionEquationOutputFile, x,step,step, false)
end
local print_grid_functions = function(self)
	vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(self.component, "lbfgs_variable")
	vtkWriter:print("q_tmp", self.q_tmp,0,0, false)
	vtkWriter:print("x_old", self.x_old,0,0, false)
	vtkWriter:print("g_old", self.g_old,0,0, false)
end

local set_gf_vectors = function(self,v_y,v_s)
	self.y=v_y
	self.s=v_s
end

local methods = {
	init=init,
	reset=reset,
	update=update,
	update_internal=update_internal,
	print_status=print_status,
	init_grid_functions=init_grid_functions,
	createGFVector=createGFVector,
	print_vector_grid_functions=print_vector_grid_functions,
	print_grid_functions=print_grid_functions,
	function_test=function_test,
	call_dot_product=call_dot_product,
	set_gf_vectors=set_gf_vectors,
	update_all_external=update_all_external,
}

local new = function()
	 local r = {
		y = {}, --stores gradient g_k
		s = {}, --stores scalar vecs. x_k
		valpha = {},
		rho = {},
		iter = 0,
		memory = 0,
		--q_tmp=GridFunction(AppSpace),
		--x_old=GridFunction(AppSpace),
		--g_old=GridFunction(AppSpace),
		q_tmp=nil,
		x_old=nil,
		g_old=nil,
		inner_gfs = {},--q_tmp,x_old,g_old
		indeces = deque.new(),
		component = nil,
		innerdomain=nil,
		boundarydomain=nil,
		dot_product=nil,
		integration_order=4,
		}
  return setmetatable(r, {__index = methods})
end

return {
  new = new,
}