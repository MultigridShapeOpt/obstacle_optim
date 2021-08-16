local deque = require "deque"

local init_grid_functions= function(self, AppSpace_1, AppSpace_2)
	if AppSpace_1 ~= nil and AppSpace_2 ~= nil then 
	print("Given Approximation Space, initialize GridFunction vectors")
		self:createGFVector(AppSpace_1, AppSpace_2)
		print("Create individual gfs")

		self.inner_gfs_1[1]=GridFunction(AppSpace_1);self.inner_gfs_1[1]:set(0.0);
		self.inner_gfs_1[2]=GridFunction(AppSpace_1);self.inner_gfs_1[2]:set(0.0);
		self.inner_gfs_1[3]=GridFunction(AppSpace_1);self.inner_gfs_1[3]:set(0.0);

		self.inner_gfs_2[1]=GridFunction(AppSpace_2);self.inner_gfs_2[1]:set(0.0);
		self.inner_gfs_2[2]=GridFunction(AppSpace_2);self.inner_gfs_2[2]:set(0.0);
		self.inner_gfs_2[3]=GridFunction(AppSpace_2);self.inner_gfs_2[3]:set(0.0);
	else 
		print("BFGS: Approximation Space is null, GridFunctions not initialized")
	end
	
end

local createGFVector = function(self, AppSpace_1, AppSpace_2)
	for i=1, self.memory do
		self.y_1[i]=GridFunction(AppSpace_1);self.y_1[i]:set(0.0);
		self.s_1[i]=GridFunction(AppSpace_1);self.s_1[i]:set(0.0);

		self.y_2[i]=GridFunction(AppSpace_2);self.y_2[i]:set(0.0);
		self.s_2[i]=GridFunction(AppSpace_2);self.s_2[i]:set(0.0);
	end
end

local reset = function(self)
	self.iter=0
	self.indices:clear()
end

local print_status = function(self)
	print("Size of y "..#self.y_1)
	print("Size of s "..#self.s_1)
	print("Value of iter "..self.iter)
	print("Value of memory "..self.memory)
	print("Size of rho "..#self.rho)
	print("Size of alpha "..#self.valpha)
end



local update = function(self,x_k_1,g_k_1,x_k_2,g_k_2)
	if self.iter>0 and self.memory>0 then

		VecScaleAssign(self.inner_gfs_1[1],1.0,g_k_1)
		VecScaleAssign(self.inner_gfs_2[1],1.0,g_k_2)
		
		VecScaleAdd2(self.inner_gfs_1[2],1.0,x_k_1,-1.0,self.inner_gfs_1[2])
		VecScaleAdd2(self.inner_gfs_2[2],1.0,x_k_2,-1.0,self.inner_gfs_2[2])
		i=self.iter%self.memory
		if i==0 then 
			i=self.memory 
		end
		
		
		VecScaleAssign(self.s_1[i],1.0,self.inner_gfs_1[2])
		VecScaleAssign(self.s_2[i],1.0,self.inner_gfs_2[2])

		self.indices:push_right(i)		
		if self.indices:length() > self.memory then 
			self.indices:pop_left()
		end
		
		VecScaleAssign(self.inner_gfs_1[2],1.0,x_k_1)
		VecScaleAssign(self.inner_gfs_2[2],1.0,x_k_2)
		
		
		VecScaleAdd2(self.inner_gfs_1[3],1.0,g_k_1,-1.0,self.inner_gfs_1[3])
		VecScaleAssign(self.y_1[i],1.0,self.inner_gfs_1[3])
		VecScaleAssign(self.inner_gfs_1[3],1.0,g_k_1)

		VecScaleAdd2(self.inner_gfs_2[3],1.0,g_k_2,-1.0,self.inner_gfs_2[3])
		VecScaleAssign(self.y_2[i],1.0,self.inner_gfs_2[3])
		VecScaleAssign(self.inner_gfs_2[3],1.0,g_k_2)

		self.rho[i]=1/self:call_dot_product(self.y_1[i],self.s_1[i], self.y_2[i],self.s_2[i])
			
		for j in self.indices:iter_right() do						
			self.valpha[j]=self.rho[j]*self:call_dot_product(self.inner_gfs_1[1],self.s_1[j], self.inner_gfs_2[1],self.s_2[j])
			VecScaleAdd2(self.inner_gfs_1[1],1.0,self.inner_gfs_1[1],-self.valpha[j],self.y_1[j])
			VecScaleAdd2(self.inner_gfs_2[1],1.0,self.inner_gfs_2[1],-self.valpha[j],self.y_2[j])
			
		end

		local gamma=self:call_dot_product(self.y_1[i],self.s_1[i], self.y_2[i],self.s_2[i])/
					self:call_dot_product(self.y_1[i],self.y_1[i], self.y_2[i],self.y_2[i])
				
		VecScaleAssign(self.inner_gfs_1[1],gamma,self.inner_gfs_1[1])
		VecScaleAssign(self.inner_gfs_2[1],gamma,self.inner_gfs_2[1])
		
		for j in self.indices:iter_left() do			
			local beta = self:call_dot_product(self.inner_gfs_1[1],self.y_1[j], self.inner_gfs_2[1],self.y_2[j])
			VecScaleAdd2(self.inner_gfs_1[1],1.0,self.inner_gfs_1[1],self.valpha[j]-beta,self.s_1[j])
			VecScaleAdd2(self.inner_gfs_2[1],1.0,self.inner_gfs_2[1],self.valpha[j]-beta,self.s_2[j])				
		end	
		VecScaleAssign(g_k_1,-1.0,self.inner_gfs_1[1])
		VecScaleAssign(g_k_2,-1.0,self.inner_gfs_2[1])
	elseif self.iter == 0 then
		print("INIT CASE")
		self.rho= {}
		VecScaleAssign(self.inner_gfs_1[2],1.0,x_k_1)
		VecScaleAssign(self.inner_gfs_1[3],1.0,g_k_1)

		VecScaleAssign(self.inner_gfs_2[2],1.0,x_k_2)
		VecScaleAssign(self.inner_gfs_2[3],1.0,g_k_2)
	else
		
	end
	self.iter=self.iter+1
end

local init = function(self,mmry,cmp_1,innerss_1,surfaces_1,scalar_product_1, cmp_2,innerss_2,surfaces_2,scalar_product_2)
	self.memory=mmry
	self.valpha[mmry]=0
	self.rho[mmry]=0
	self.component_1=cmp_1 
	self.innerdomain_1=innerss_1
	self.boundarydomain_1=surfaces_1
	self.dot_product_1=scalar_product_1
	self.component_2=cmp_2 
	self.innerdomain_2=innerss_2
	self.boundarydomain_2=surfaces_2
	self.dot_product_2=scalar_product_2
end

local call_dot_product= function(self,gf_a_1,gf_b_1, gf_a_2,gf_b_2)
		return self.dot_product_1(gf_a_1,gf_b_1,self.component_1,self.boundarydomain_1,self.innerdomain_1,self.integration_order) 
			 + self.dot_product_2(gf_a_2,gf_b_2,self.component_2,self.innerdomain_2,self.integration_order) 
end

local methods = {
	init=init,
	reset=reset,
	print_status=print_status,
	init_grid_functions=init_grid_functions,
	createGFVector=createGFVector,
	call_dot_product=call_dot_product,
	update=update,

}

local new = function()
	 local r = {
		y_1 = {}, 
		s_1 = {}, 
		y_2 = {}, 
		s_2 = {}, 
		valpha = {},
		rho = {},
		iter = 0,
		memory = 0,

		inner_gfs_1 = {},
		indices = deque.new(),
		component_1 = nil,
		innerdomain_1=nil,
		boundarydomain_1=nil,
		dot_product_1=nil,
		integration_order=4,

		inner_gfs_2 = {},
		component_2 = nil,
		innerdomain_2=nil,
		boundarydomain_2=nil,
		dot_product_2=nil,
		}
  return setmetatable(r, {__index = methods})
end

return {
  new = new,
}